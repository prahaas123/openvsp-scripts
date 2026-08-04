import os
import csv
import glob
import uuid
import math
import numpy as np
import openvsp as vsp  # type: ignore
from scipy.optimize import differential_evolution, NonlinearConstraint
import pyvista as pv  # type: ignore

VELOCITY = 12.0          # [m/s] cruise airspeed
RHO = 1.225
MU = 1.81e-5
A_SOUND = 343.0          # [m/s]
Q_INF = 0.5 * RHO * VELOCITY ** 2
 
AIRFOIL_FILE = os.path.join("Airfoils", "mh61.dat")

# Design Variable Bounds
VAR_NAMES = ["root_chord", "taper", "sweep", "washout", "span", "dihedral", "alpha"]
BOUNDS = [
    (0.15, 0.40),    # root_chord [m]
    (0.25, 1.00),    # taper ratio
    (0.0, 45.0),     # LE sweep [deg]
    (-4.0, 0.0),     # washout [deg]
    (0.50, 1.10),    # full span [m]
    (0.0, 5.0),      # dihedral [deg]
    (2.0, 8.0),      # alpha [deg]
]

# Weight Estimate
M_FIXED = 0.265          # [kg]      avionics, battery, propulsion, payload
M_SPAR_PER_M = 0.044     # [kg/m]    spar mass per unit span
M_SKIN_PER_M2 = 0.950    # [kg/m^2]  3D-printed skin + ribs per unit planform area

def total_weight(span, s_ref):
    return 9.81 * (M_FIXED + M_SPAR_PER_M * span + M_SKIN_PER_M2 * s_ref)

# Constraint Targets
STATIC_MARGIN = 0.05
CM_TOL = 0.02
AR_MIN = 2.0
AR_MAX = 6.0
TIP_CHORD_MIN = 0.05

ALPHA_STEP_DEG = 1.0

# Pentalty Weights
PENALTY_WEIGHT = 20.0
FAILED_RUN_OBJECTIVE = 1.0e6

# Parasitic Drag Model
LAMINAR_FRACTION = 0.40  # fraction of chord assumed laminar
FF_MACH_MIN = 0.30       # below this Mach the compressibility form-factor term is skipped
TC_DEFAULT = 0.10        # fallback t/c
XTC_DEFAULT = 0.30       # fallback chordwise location of max t/c

# VSP Settings
WING_SPAN_RES = 10
WING_CHORD_RES = 25
WAKE_NUM_ITER = 5
NCPU = 8

# Logging
LOG_CSV = "optimization_results.csv"
LOG_FIELDS = [
    "run_id", "objective",
    "root_chord", "taper", "sweep", "washout", "span", "dihedral", "alpha",
    "AR", "mass_kg", "weight", "lift", "lift_margin",
    "CL", "CD", "LD", "drag_N",
    "CM_cg", "x_cg", "SM",
]

def init_log():
    with open(LOG_CSV, "w", newline="") as f:
        csv.DictWriter(f, fieldnames=LOG_FIELDS).writeheader()

def append_log(row):
    with open(LOG_CSV, "a", newline="") as f:
        csv.DictWriter(f, fieldnames=LOG_FIELDS, extrasaction="ignore").writerow(row)

def lookup_best():
    best = None
    if not os.path.isfile(LOG_CSV):
        return None
    with open(LOG_CSV, newline="") as f:
        for row in csv.DictReader(f):
            try:
                obj = float(row["objective"])
            except (KeyError, ValueError):
                continue
            if best is None or obj < float(best["objective"]):
                best = row
    return best

# Geometric Evaluations
def planform(root_chord, taper, span):
    tip_chord = root_chord * taper
    s_ref = 0.5 * (root_chord + tip_chord) * span
    ar = span ** 2 / max(s_ref, 1e-12)
    mac = (2.0 / 3.0) * root_chord * (1.0 + taper + taper ** 2) / (1.0 + taper)
    y_mac = (span / 2.0) * (1.0 + 2.0 * taper) / (3.0 * (1.0 + taper))
    return dict(tip_chord=tip_chord, s_ref=s_ref, ar=ar, mac=mac, y_mac=y_mac)

def x_ac_estimate(root_chord, taper, span, sweep_deg):
    g = planform(root_chord, taper, span)
    return g["y_mac"] * math.tan(math.radians(sweep_deg)) + 0.25 * g["mac"]

def airfoil_thickness(path):
    try:
        pts = []
        with open(path, "r") as f:
            for line in f:
                parts = line.replace(",", " ").split()
                if len(parts) < 2:
                    continue
                try:
                    x, y = float(parts[0]), float(parts[1])
                except ValueError:
                    continue
                if -0.05 <= x <= 1.05 and abs(y) <= 0.6:
                    pts.append((x, y))
        if len(pts) < 20:
            raise ValueError("too few coordinate pairs")
 
        pts = np.array(pts)
        i_le = int(np.argmin(pts[:, 0]))
        upper = pts[: i_le + 1][::-1]
        lower = pts[i_le:]
        if len(upper) < 5 or len(lower) < 5:
            raise ValueError("could not split surfaces")
        if np.mean(upper[:, 1]) < np.mean(lower[:, 1]):
            upper, lower = lower, upper
 
        xq = np.linspace(0.0, 1.0, 200)
        yu = np.interp(xq, np.sort(upper[:, 0]), upper[np.argsort(upper[:, 0]), 1])
        yl = np.interp(xq, np.sort(lower[:, 0]), lower[np.argsort(lower[:, 0]), 1])
        thick = yu - yl
        i_max = int(np.argmax(thick))
        tc = float(thick[i_max])
        xtc = float(xq[i_max])
        if not (0.02 < tc < 0.4):
            raise ValueError(f"implausible t/c = {tc}")
        return tc, xtc
    except Exception as e:
        print(f"[warn] airfoil thickness parse failed ({e}); using t/c={TC_DEFAULT}")
        return TC_DEFAULT, XTC_DEFAULT

TC_MAX, X_TC_MAX = airfoil_thickness(AIRFOIL_FILE)

def estimate_cd0(mac, s_ref, sweep_deg):
    re = RHO * VELOCITY * mac / MU
    re = max(re, 1.0e4)
    mach = VELOCITY / A_SOUND
 
    cf_lam = 1.328 / math.sqrt(re)
    cf_turb = 0.455 / (math.log10(re) ** 2.58) / ((1.0 + 0.144 * mach ** 2) ** 0.65)
    cf = LAMINAR_FRACTION * cf_lam + (1.0 - LAMINAR_FRACTION) * cf_turb
 
    # Raymer component form factor for a lifting surface
    ff = (1.0 + (0.6 / max(X_TC_MAX, 0.05)) * TC_MAX + 100.0 * TC_MAX ** 4)
    if mach >= FF_MACH_MIN:
        ff *= max(1.0, 1.34 * (mach ** 0.18) * (math.cos(math.radians(sweep_deg)) ** 0.28))
 
    s_wet = 2.0 * s_ref * (1.0 + 1.2 * TC_MAX)
    return float(cf * ff * s_wet / max(s_ref, 1e-12))

# VSP Helper Functions
_ANALYSIS_INPUTS = {}
CDO_INPUT_NAME = None

def analysis_inputs(analysis_name):
    if analysis_name not in _ANALYSIS_INPUTS:
        try:
            _ANALYSIS_INPUTS[analysis_name] = list(vsp.GetAnalysisInputNames(analysis_name))
        except Exception:
            _ANALYSIS_INPUTS[analysis_name] = []
    return _ANALYSIS_INPUTS[analysis_name]

def find_input(analysis_name, candidates):
    available = analysis_inputs(analysis_name)
    lowered = {n.lower(): n for n in available}
    for c in candidates:
        if c.lower() in lowered:
            return lowered[c.lower()]
    return None

def set_double(analysis, candidates, value):
    name = find_input(analysis, candidates)
    if name is not None:
        vsp.SetDoubleAnalysisInput(analysis, name, [float(value)])
    return name

def set_int(analysis, candidates, value):
    name = find_input(analysis, candidates)
    if name is not None:
        vsp.SetIntAnalysisInput(analysis, name, [int(value)])
    return name

# Geometry Generation
def generate_wing(name, span, root_chord, taper, sweep_deg, dihedral_deg, washout_deg, airfoil_file=AIRFOIL_FILE, export_stl=False):
    tip_chord = root_chord * taper
 
    vsp.VSPCheckSetup()
    vsp.ClearVSPModel()
    wing_id = vsp.AddGeom("WING")
    vsp.SetParmVal(wing_id, "TotalSpan",      "WingGeom", span)
    vsp.SetParmVal(wing_id, "Root_Chord",     "XSec_1",   root_chord)
    vsp.SetParmVal(wing_id, "Tip_Chord",      "XSec_1",   tip_chord)
    vsp.SetParmVal(wing_id, "Sweep",          "XSec_1",   sweep_deg)
    vsp.SetParmVal(wing_id, "Dihedral",       "XSec_1",   dihedral_deg)
    vsp.SetParmVal(wing_id, "Twist",          "XSec_1",   washout_deg)
    vsp.SetParmVal(wing_id, "Twist_Location", "XSec_1",   1.0)
    vsp.SetParmVal(wing_id, "SectTess_U",     "XSec_1",   WING_SPAN_RES)
    vsp.SetParmVal(wing_id, "Tess_W",         "Shape",    WING_CHORD_RES)
 
    # Airfoil selection
    root_xsec_surf = vsp.GetXSecSurf(wing_id, 0)
    vsp.ChangeXSecShape(root_xsec_surf, 0, vsp.XS_FILE_AIRFOIL)
    root_xsec = vsp.GetXSec(root_xsec_surf, 0)
    vsp.ReadFileAirfoil(root_xsec, airfoil_file)
 
    tip_xsec_surf = vsp.GetXSecSurf(wing_id, 1)
    vsp.ChangeXSecShape(tip_xsec_surf, 1, vsp.XS_FILE_AIRFOIL)
    tip_xsec = vsp.GetXSec(tip_xsec_surf, 1)
    vsp.ReadFileAirfoil(tip_xsec, airfoil_file)
 
    vsp.SetSetFlag(wing_id, 1, True)
    vsp.Update()
 
    vsp3_path = f"{name}.vsp3"
    vsp.WriteVSPFile(vsp3_path)
 
    stl_path = None
    if export_stl:
        stl_path = f"{name}.stl"
        vsp.ExportFile(stl_path, 0, vsp.EXPORT_STL)
 
    return vsp3_path, stl_path

# VSPAero Runs
def _result_array(res_id, names):
    for n in names:
        try:
            vals = vsp.GetDoubleResults(res_id, n)
            if vals is not None and len(vals) > 0:
                return np.asarray(vals, dtype=float)
        except Exception:
            continue
    raise KeyError(f"None of {names} found in VSPAERO results.")

def _run_vspaero(vsp3_path, s_ref, b_ref, c_ref, x_cg, cd0,
                 alpha_start, alpha_end, alpha_npts,
                 beta_start=0.0, beta_end=0.0, beta_npts=1):
    mach = VELOCITY / A_SOUND
 
    vsp.ClearVSPModel()
    vsp.ReadVSPFile(vsp3_path)
 
    geom_analysis = "VSPAEROComputeGeometry"
    vsp.SetAnalysisInputDefaults(geom_analysis)
    set_int(geom_analysis, ["GeomSet"], vsp.SET_NONE)
    set_int(geom_analysis, ["ThinGeomSet"], vsp.SET_ALL)
    vsp.ExecAnalysis(geom_analysis)
 
    aero = "VSPAEROSweep"
    vsp.SetAnalysisInputDefaults(aero)
 
    set_double(aero, ["Sref"], s_ref)
    set_double(aero, ["bref"], b_ref)
    set_double(aero, ["cref"], c_ref)
    set_double(aero, ["Xcg"], x_cg)
    set_double(aero, ["Ycg"], 0.0)
    set_double(aero, ["Zcg"], 0.0)
 
    set_double(aero, ["AlphaStart"], alpha_start)
    set_double(aero, ["AlphaEnd"], alpha_end)
    set_int(aero, ["AlphaNpts"], alpha_npts)
 
    set_double(aero, ["BetaStart"], beta_start)
    set_double(aero, ["BetaEnd"], beta_end)
    set_int(aero, ["BetaNpts"], beta_npts)
 
    set_double(aero, ["MachStart"], mach)
    set_double(aero, ["MachEnd"], mach)
    set_int(aero, ["MachNpts"], 1)
    set_double(aero, ["Vinf"], VELOCITY)
    set_double(aero, ["Rho"], RHO)
 
    set_int(aero, ["Symmetry"], 0)
 
    set_int(aero, ["WakeNumIter"], WAKE_NUM_ITER)
    set_int(aero, ["NCPU"], NCPU)
 
    if CDO_INPUT_NAME is not None:
        vsp.SetDoubleAnalysisInput(aero, CDO_INPUT_NAME, [float(cd0)])
 
    name = find_input(aero, ["RedirectFile"])
    if name is not None:
        vsp.SetStringAnalysisInput(aero, name, [f"{vsp3_path}_log.txt"])
 
    vsp.ExecAnalysis(aero)
 
    res_id = vsp.FindLatestResultsID("VSPAERO_Polar")
    out = {
        "CL": _result_array(res_id, ["CLtot", "CL"]),
        "CD": _result_array(res_id, ["CDtot", "CD"]),
        "CMy": _result_array(res_id, ["CMytot", "CMy", "CMmtot"]),
    }
    try:
        out["CMx"] = _result_array(res_id, ["CMxtot", "CMx", "CMltot"])
        out["CMz"] = _result_array(res_id, ["CMztot", "CMz", "CMntot"])
    except KeyError:
        out["CMx"] = np.zeros_like(out["CL"])
        out["CMz"] = np.zeros_like(out["CL"])
    try:
        out["Alpha"] = _result_array(res_id, ["Alpha", "AoA"])
    except KeyError:
        out["Alpha"] = np.linspace(alpha_start, alpha_end, alpha_npts)
    try:
        out["Beta"] = _result_array(res_id, ["Beta"])
    except KeyError:
        out["Beta"] = np.linspace(beta_start, beta_end, beta_npts)
 
    if CDO_INPUT_NAME is None:
        out["CD"] = out["CD"] + float(cd0)
 
    return out

def _pick(values, coord, target):
    coord = np.asarray(coord, dtype=float)
    n = min(len(coord), len(values))
    if n == 0:
        return float("nan")
    return float(np.asarray(values, dtype=float)[: n][int(np.argmin(np.abs(coord[:n] - target)))])

# Objective Evaluation
def evaluate_geometry(x):
    root_chord, taper, sweep, washout, span, dihedral, alpha = x
    g = planform(root_chord, taper, span)
    return np.array([g["ar"], g["tip_chord"]])

def _quad_pen(violation, scale):
    if violation <= 0.0:
        return 0.0
    return PENALTY_WEIGHT * (violation / max(abs(scale), 1e-9)) ** 2

def evaluate_aero_objective(x):
    root_chord, taper, sweep, washout, span, dihedral, alpha = [float(v) for v in x]
    run_id = f"wing_{uuid.uuid4().hex[:8]}"
 
    try:
        g = planform(root_chord, taper, span)
        s_ref, ar, mac = g["s_ref"], g["ar"], g["mac"]
        weight = total_weight(span, s_ref)
        cd0 = estimate_cd0(mac, s_ref, sweep)
 
        vsp3_path, _ = generate_wing(run_id, span, root_chord, taper,
                                     sweep, dihedral, washout)
 
        # Alpha sweep, moment reference at the root LE (x=0)
        r1 = _run_vspaero(
            vsp3_path, s_ref, span, mac, x_cg=0.0, cd0=cd0,
            alpha_start=alpha - ALPHA_STEP_DEG,
            alpha_end=alpha + ALPHA_STEP_DEG,
            alpha_npts=3,
        )
 
        a = np.asarray(r1["Alpha"], dtype=float)
        cl_m = _pick(r1["CL"], a, alpha - ALPHA_STEP_DEG)
        cl_0 = _pick(r1["CL"], a, alpha)
        cl_p = _pick(r1["CL"], a, alpha + ALPHA_STEP_DEG)
        cm_m = _pick(r1["CMy"], a, alpha - ALPHA_STEP_DEG)
        cm_0 = _pick(r1["CMy"], a, alpha)
        cm_p = _pick(r1["CMy"], a, alpha + ALPHA_STEP_DEG)
        cd_0 = _pick(r1["CD"], a, alpha)
 
        d_alpha_rad = 2.0 * math.radians(ALPHA_STEP_DEG)
        cl_alpha = (cl_p - cl_m) / d_alpha_rad
        cm_alpha = (cm_p - cm_m) / d_alpha_rad
 
        if not np.isfinite(cl_alpha) or abs(cl_alpha) < 1e-6:
            raise RuntimeError("degenerate CL_alpha")
 
        # Neutral point about the root LE, then CG at the target static margin.
        x_np = -mac * (cm_alpha / cl_alpha)
        x_cg = x_np - STATIC_MARGIN * mac
 
        x_np_est = x_ac_estimate(root_chord, taper, span, sweep)
        if not np.isfinite(x_np) or abs(x_np - x_np_est) > 0.5 * mac:
            print(f"  [{run_id}] [warn] NP {x_np:.4f} m vs quarter-MAC estimate "
                  f"{x_np_est:.4f} m -- check CMy sign convention / cref")
 
        # Moment about the CG at the design alpha (shift theorem, no extra run).
        cm_cg = cm_0 + cl_0 * x_cg / mac
        sm_actual = (x_np - x_cg) / mac
 
        # Forces
        lift = Q_INF * s_ref * cl_0
        drag = Q_INF * s_ref * cd_0
        ld = cl_0 / cd_0 if cd_0 > 1e-9 else float("nan")
 
        # Penalties
        pen_lift = _quad_pen(weight - lift, weight)
        pen_trim = _quad_pen(abs(cm_cg) - CM_TOL, max(CM_TOL, 1e-3))
        penalty = pen_lift + pen_trim
        feasible = penalty <= 0.0
        objective = drag + penalty
        
        if feasible:
            append_log({
                "run_id": run_id,
                "objective": round(objective, 6), "taper": round(taper, 6),
                "sweep": round(sweep, 6), "washout": round(washout, 6),
                "span": round(span, 6), "dihedral": round(dihedral, 6),
                "alpha": round(alpha, 6), "AR": round(ar, 6),
                "tip_chord": round(g["tip_chord"], 6),
                "weight": round(weight, 6), "mass_kg": round(weight/9.81, 6),
                "lift": round(lift, 6), "lift_margin": round(lift - weight, 6),
                "CL": round(cl_0, 6), "CD": round(cd_0, 6),
                "LD": round(ld, 6), "drag_N": round(drag, 6),
                "CM_cg": round(cm_cg, 6),
                "x_cg": round(x_cg, 6), "SM": round(sm_actual, 6),
            })
 
        tag = "FEASIBLE" if feasible else "infeasible"
        print(f"  [{run_id}] {tag:10s} D={drag:.4f} N  obj={objective:.4f}  "
              f"L/D={ld:.2f}  L-W={lift - weight:+.3f} N  CMcg={cm_cg:+.4f}  "
              f"m={weight/9.81:.3f} kg  AR={ar:.2f}")
 
        return float(objective)
 
    except Exception as e:
        print(f"  [{run_id}] FAILED: {e}")
        return FAILED_RUN_OBJECTIVE
 
    finally:
        for filename in glob.glob(f"{run_id}*"):
            try:
                os.remove(filename)
            except OSError:
                pass

def visualize_stl(stl_path):
    if pv is None:
        print("[info] pyvista not available; skipping visualization.")
        return
    if not os.path.exists(stl_path):
        print("Error: STL not found.")
        return
    mesh = pv.read(stl_path)
    plotter = pv.Plotter(title="Optimized Wing")
    plotter.add_mesh(mesh, color="lightblue", show_edges=True, smooth_shading=True)
    plotter.add_axes()
    plotter.add_floor(face="-z", i_resolution=10, j_resolution=10,
                      color="gray", opacity=0.2)
    print("Opening PyVista window...")
    plotter.show()

def main():
    init_log()
 
    print(f"[setup] airfoil t/c = {TC_MAX:.4f} at x/c = {X_TC_MAX:.3f}")
    print(f"[setup] q_inf = {Q_INF:.2f} Pa at V = {VELOCITY:.1f} m/s")
    print("Starting SciPy DE optimization...\n")
 
    geom_constraint = NonlinearConstraint(
        evaluate_geometry,
        lb=[AR_MIN, TIP_CHORD_MIN],
        ub=[AR_MAX, np.inf],
    )
 
    result = None
    try:
        result = differential_evolution(
            evaluate_aero_objective,
            bounds=BOUNDS,
            constraints=(geom_constraint,),
            strategy="rand1bin",
            recombination=0.9,
            mutation=(0.5, 1.0),
            popsize=10,
            maxiter=30,
            tol=1e-3,
            init="latinhypercube",
            polish=False,
            seed=1,
            disp=True,
        )
        print("\nOptimization finished.")
    except KeyboardInterrupt:
        print("\n\nInterrupted. Reporting the best feasible logged design.")
 
    row = lookup_best()
    if row is None:
        print("No feasible designs were logged.")
        if result is not None:
            print(f"DE returned x = {result.x}, f = {result.fun}")
        return
 
    p = {k: float(row[k]) for k in
         ["root_chord", "taper", "sweep", "washout", "span", "dihedral", "alpha"]}
 
    print("\n--- OPTIMAL WING ---")
    print("\n--- PARAMETERS ---")
    for k, v in p.items():
        print(f"{k:12s}: {v:.4f}")
 
    print("\n--- GEOMETRY ---")
    for k in ["S_ref", "AR", "MAC", "tip_chord"]:
        print(f"{k:12s}: {float(row[k]):.4f}")
    print(f"{'x_np':12s}: {float(row['x_np']):.4f} m aft of root LE")
    print(f"{'x_cg':12s}: {float(row['x_cg']):.4f} m  (SM = {float(row['SM_actual']) * 100:.2f}% MAC)")
 
    print("\n--- AERODYNAMICS ---")
    for k in ["CL", "CD", "CD0", "LD", "drag_N", "CM_cg"]:
        print(f"{k:12s}: {float(row[k]):.5f}")
 
    print("\n--- WEIGHT ---")
    print(f"{'mass':12s}: {float(row['mass_kg']):.4f} kg")
    print(f"{'weight':12s}: {float(row['weight']):.4f} N")
    print(f"{'lift':12s}: {float(row['lift']):.4f} N")
    print(f"{'margin':12s}: {float(row['lift_margin']):+.4f} N")
 
    print(f"\nAll designs logged to: {LOG_CSV}\n")
 
    _, stl_path = generate_wing(
        "Optimized_Wing", p["span"], p["root_chord"], p["taper"],
        p["sweep"], p["dihedral"], p["washout"], export_stl=True,
    )
    if stl_path:
        visualize_stl(stl_path)

if __name__ == "__main__":
    main()