import os
import openvsp as vsp # type: ignore
import csv
import pyvista as pv
import numpy as np
import glob
import uuid
from scipy.optimize import differential_evolution, NonlinearConstraint

vsp_exe = r"C:\Program Files\OpenVSP-3.47.0\vsp.exe"

wing_span_res = 10
wing_chord_res = 25
velocity = 20 # m/s
alpha = 3 # degrees AoA

airfoil_file = r"Airfoils\pw75.dat"

MAX_WEIGHT = 9.81 * 0.5   # Newtons
MIN_WEIGHT = 9.81 * 0.3   # Newtons
WING_LOADING = 40  # N/m^2
STATIC_MARGIN = 0.05
CM_MIN = -0.05   # lower bound on CM about CG
CM_MAX =  0.05   # upper bound on CM about CG
AR_MIN = 2.0
AR_MAX = 6.0
TIP_CHORD_MIN = 0.05
MIN_S_REF = MIN_WEIGHT / WING_LOADING  # m2
MAX_S_REF = MAX_WEIGHT / WING_LOADING  # m2

LOG_CSV = "optimization_results.csv"
LOG_FIELDS = [
    "run_id",
    "root_chord", "taper", "sweep", "twist", "span",
    "LD", "CL", "CD", "CM_cg",
    "AR", "lift", "x_cg"
]

def main():   
    init_log()
    print("Starting SciPy DE Optimization...")
    bounds = [(0.18, 0.3),   # Root chord
              (0.6, 0.95),   # Taper ratio
              (20.0, 40.0), # Sweep angle
              (-5.0, 0.0), # Washout angle
              (0.5, 0.65)]   # Wingspan

    # Early-rejection geometric constraints
    geom_constraint = NonlinearConstraint(
        evaluate_geometry, 
        lb=[MIN_S_REF, AR_MIN, TIP_CHORD_MIN], 
        ub=[MAX_S_REF, AR_MAX, np.inf]
    )

    # Run SciPy DE
    try:
        result = differential_evolution(
            evaluate_aero_objective,
            bounds=bounds,
            constraints=(geom_constraint,),
            strategy='rand1bin',
            recombination=0.9,
            mutation=(0.5, 1.0),
            popsize=6,               # Population members = popsize * parameters 
            maxiter=30,              # n_max_gen
            tol=1e-3,                # ftol
            seed=1,
            disp=True
        )
        print(f"\nOptimization finished.")
    except KeyboardInterrupt:
        print("\n\nOptimization interrupted by user! Generating results from the best logged design so far...")

    # Look up values from the log (safest way to get the best strictly feasible run)
    row = lookup_best()
    if row is None:
        print("No feasible designs were logged — cannot print results summary.")
        return
        
    best_root   = float(row["root_chord"])
    best_taper  = float(row["taper"])
    best_sweep  = float(row["sweep"])
    best_twist  = float(row["twist"])
    best_span   = float(row["span"])
    ld          = float(row["LD"])
    cm          = float(row["CM_cg"])
    ar          = float(row["AR"])
    lift        = float(row["lift"])
    x_cg        = float(row["x_cg"])

    print("\n--- OPTIMAL WING GEOMETRY ---")
    print(f"\n--- AERODYNAMICS ---")
    print(f"L/D         : {ld:.4f}")
    print(f"CL          : {float(row['CL']):.4f}")
    print(f"CD          : {float(row['CD']):.4f}")
    print(f"CM_cg       : {cm:.4f}")
    print(f"Total lift  : {lift:.4f} N")
    print(f"\n--- GEOMETRY ---")
    print(f"Aspect ratio: {ar:.4f}  (min: {AR_MIN}  max: {AR_MAX})")
    print(f"CG          : {x_cg:.4f} m  (aft of root LE,  SM={STATIC_MARGIN*100:.0f}% MAC)")
    print(f"\n--- PARAMETERS ---")
    print(f"Params      : Root={best_root:.4f}  Taper={best_taper:.4f}  Sweep={best_sweep:.4f}  Twist={best_twist:.4f}  Span={best_span:.4f}")
    print(f"\nFeasible designs logged to: {LOG_CSV}\n")
    
    # Generate STL
    stl_path, _ = generate_wing("Optimized_Wing", best_span * 1000, best_root * 1000, best_taper, best_sweep, 0.0, best_twist, airfoil_file)
    visualize_stl(stl_path)
    
def evaluate_geometry(x):
    root_chord, taper, sweep, twist, span = x
    Sref = 0.5 * (root_chord + root_chord * taper) * span
    AR = aspect_ratio(root_chord, taper, span)
    tip_chord = root_chord * taper
    return np.array([Sref, AR, tip_chord])

def evaluate_aero_objective(x):
    root_chord, taper, sweep, twist, span = x
    run_id = f"wing_{uuid.uuid4().hex[:8]}" 
    
    try:
        stl_path, analysis_path = generate_wing(run_id, span, root_chord, taper, sweep, 0.0, twist, airfoil_file)
        Sref = 0.5 * (root_chord + root_chord * taper) * span
        x_cg = calc_cg(root_chord, taper, span, sweep)
        CL, CD, LD, CM_cg  = vsp_point(analysis_path, velocity, alpha, Sref, span, root_chord, x_cg)
        lift  = 0.5 * CL * 1.225 * velocity**2 * Sref
        AR    = aspect_ratio(root_chord, taper, span)
        
        # Check Aerodynamic Constraints
        penalty = 0
        is_feasible = True
        if lift < MAX_WEIGHT:
            penalty += 1e5 * (MAX_WEIGHT - lift) # Scale penalty by how badly it failed
            is_feasible = False
        if CM_cg < CM_MIN or CM_cg > CM_MAX:
            penalty += 1e5
            is_feasible = False

        if is_feasible:
            print(f"  [{run_id}] FEASIBLE! L/D={LD:.3f}")
            append_log({
                "run_id"      : run_id,
                "root_chord"  : round(float(root_chord),  6),
                "taper"       : round(float(taper),       6),
                "sweep"       : round(float(sweep),       6),
                "twist"       : round(float(twist),       6),
                "span"        : round(float(span),        6),
                "LD"          : round(LD,                 6),
                "CL"          : round(CL,                 6),
                "CD"          : round(CD,                 6),
                "CM_cg"       : round(CM_cg,              6),
                "AR"          : round(AR,                 6),
                "lift"        : round(lift,               6),
                "x_cg"        : round(x_cg,               6),
            })
        else:
            print(f"  [{run_id}] Failed Aero Constraints. Lift={lift:.2f}N, CM={CM_cg:.3f}")

        return -(CL ** 2 / CD) + penalty

    except Exception as e:
        print(f"Run {run_id} failed: {e}")
        return 1e10 # Massive penalty for crashed runs
    finally:
        for filename in glob.glob(f"{run_id}*"):
            try:
                os.remove(filename)
            except OSError:
                pass

def generate_wing(wing_name, wingspan, root_chord, taper_ratio, sweep_angle, dihedral_angle, twist_angle, airfoil_file):
    tip_chord = root_chord * taper_ratio
    airfoil_fwd = airfoil_file.replace("\\", "/")

    vsp.VSPCheckSetup()
    vsp.ClearVSPModel()
    wing_id = vsp.AddGeom( "WING" )
    vsp.SetParmVal( wing_id, "TotalSpan",      "WingGeom", wingspan )
    vsp.SetParmVal( wing_id, "Root_Chord",     "XSec_1",   root_chord )
    vsp.SetParmVal( wing_id, "Tip_Chord",      "XSec_1",   tip_chord )
    vsp.SetParmVal( wing_id, "Sweep",          "XSec_1",   sweep_angle )
    vsp.SetParmVal( wing_id, "Dihedral",       "XSec_1",   dihedral_angle )
    vsp.SetParmVal( wing_id, "Twist",          "XSec_1",   twist_angle )
    vsp.SetParmVal( wing_id, "Twist_Location", "XSec_1",   1.0 )
    vsp.SetParmVal( wing_id, "SectTess_U",     "XSec_1",   wing_span_res )
    vsp.SetParmVal( wing_id, "Tess_W",         "Shape",    wing_chord_res )
    
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

    # Finalize and export
    vsp.Update()
    stl_path = f"{wing_name}.stl"
    analysis_path = f"{wing_name}.vsp3"
    vsp.WriteVSPFile(analysis_path)
    vsp.ExportFile(stl_path, 0, vsp.EXPORT_STL)

    return stl_path, analysis_path

def visualize_stl(stl_path):
    if os.path.exists(stl_path):
        mesh = pv.read(stl_path)
        plotter = pv.Plotter(title="WatArrow Delta Wing")
        plotter.add_mesh(mesh, color="lightblue", show_edges=True, smooth_shading=True)
        plotter.add_axes()
        plotter.add_floor(face='-z', i_resolution=10, j_resolution=10, color='gray', opacity=0.2)
        print("Opening PyVista window...")
        plotter.show()
    else:
        print("Error: STL not found.")

def vsp_point(vsp3_path, vin, alpha, Sref, bref, cref, x_cg):
    mach = vin / 343.0

    # Meshing    
    vsp.ClearVSPModel()
    vsp.ReadVSPFile(vsp3_path)
    geom_analysis = "VSPAEROComputeGeometry"
    vsp.SetAnalysisInputDefaults(geom_analysis)
    vsp.SetIntAnalysisInput(geom_analysis, "GeomSet", [vsp.SET_NONE])      
    vsp.SetIntAnalysisInput(geom_analysis, "ThinGeomSet", [vsp.SET_ALL])
    vsp.ExecAnalysis(geom_analysis)
    
    # Aero Analysis
    aero_analysis = "VSPAEROSweep"
    vsp.SetAnalysisInputDefaults(aero_analysis)
    vsp.SetDoubleAnalysisInput(aero_analysis, "Sref", [Sref])
    vsp.SetDoubleAnalysisInput(aero_analysis, "cref", [cref])
    vsp.SetDoubleAnalysisInput(aero_analysis, "bref", [bref])
    vsp.SetDoubleAnalysisInput(aero_analysis, "AlphaStart", [alpha])
    vsp.SetIntAnalysisInput(aero_analysis, "AlphaNpts", [1])
    vsp.SetDoubleAnalysisInput(aero_analysis, "AlphaEnd", [alpha])
    vsp.SetDoubleAnalysisInput(aero_analysis, "MachStart", [mach])
    vsp.SetIntAnalysisInput(aero_analysis, "MachNpts", [1])
    vsp.SetIntAnalysisInput(aero_analysis, "WakeNumIter", [6]) 
    vsp.SetDoubleAnalysisInput(aero_analysis, "Vinf", [vin])
    vsp.SetDoubleAnalysisInput(aero_analysis, "Xcg", [x_cg])
    vsp.SetIntAnalysisInput(aero_analysis, "NCPU", [8])
    vsp.SetStringAnalysisInput(aero_analysis, "RedirectFile", [f"{vsp3_path}_log.txt"])
    rid = vsp.ExecAnalysis(aero_analysis)

    # Results
    polar_res = vsp.FindLatestResultsID("VSPAERO_Polar")
    cl = vsp.GetDoubleResults(polar_res, "CLtot")[0]
    cd = vsp.GetDoubleResults(polar_res, "CDtot")[0]
    cm = vsp.GetDoubleResults(polar_res, "CMytot")[0]
    return cl, cd, cl/cd, cm

def estimate_wetted_area(root_chord, taper, span):
    tip_chord     = root_chord * taper
    planform_area = 0.5 * (root_chord + tip_chord) * span
    return 2.0 * planform_area * 1.02

def calc_cg(root_chord, taper, span, sweep_angle, static_margin=STATIC_MARGIN):
    mac   = (2.0 / 3.0) * root_chord * (1 + taper + taper**2) / (1 + taper)
    y_mac = (span / 2.0) * (1 + 2 * taper) / (3 * (1 + taper))
    x_ac  = y_mac * np.tan(np.radians(sweep_angle)) + 0.25 * mac
    x_cg  = x_ac - static_margin * mac
    return x_cg

def aspect_ratio(root_chord, taper, span):
    tip_chord     = root_chord * taper
    planform_area = 0.5 * (root_chord + tip_chord) * span
    return span**2 / planform_area

def stall_speed(cl, sref, w=MAX_WEIGHT, rho=1.225):
    return np.sqrt((2*w) / (rho*sref*cl))

def root_bending_moment(root_chord, taper, span, CL, vin=velocity, rho=1.225):
    q = 0.5 * rho * vin**2
    s = span / 2.0
    return q * CL * root_chord * s**2 * (1 + 2 * taper) / 6

def init_log():
    with open(LOG_CSV, "w", newline="") as f:
        csv.DictWriter(f, fieldnames=LOG_FIELDS).writeheader()

def append_log(row: dict):
    with open(LOG_CSV, "a", newline="") as f:
        csv.DictWriter(f, fieldnames=LOG_FIELDS, extrasaction="ignore").writerow(row)

def lookup_best():
    best = None
    with open(LOG_CSV, newline="") as f:
        for row in csv.DictReader(f):
            if best is None or float(row["LD"]) > float(best["LD"]):
                best = row
    return best

# AngelScript array helpers
def iarr(v):
    return f"array<int> = {{{v}}}"

def darr(v):
    return f"array<double> = {{{v}}}"

if __name__ == "__main__":
    main()