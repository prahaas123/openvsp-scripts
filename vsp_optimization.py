import os
import csv
import subprocess
import pyvista as pv
import numpy as np
import glob
import uuid
from pymoo.algorithms.soo.nonconvex.de import DE
from pymoo.core.problem import ElementwiseProblem
from pymoo.optimize import minimize
from pymoo.termination.default import DefaultSingleObjectiveTermination

vsp_exe = r"C:\Program Files\OpenVSP-3.47.0\vsp.exe"

wing_span_res = 20
wing_chord_res = 50
velocity = 20 # m/s
alpha = 5 # degrees AoA

airfoil_file = r"C:\Users\kprah\Desktop\Prahaas\WatArrow\CFD Automation\Airfoils\dae21.dat"

SCORING = {
    # term_name     : (weight,  reference_value)
    "ld_ratio"      : (0.45,     15.0),  # typical/target L/D
    "wetted_area"   : (0.2,      0.1),   # reference wetted area (m²)
    "stall_speed"   : (0.25,    10.0),   # reference stall speed (m/s)
    "bending_moment": (0.1,     1.5),    # reference RBM (N·m)
}

MAX_WEIGHT = 6   # Newtons
STATIC_MARGIN = 0.05
CM_MIN = -0.08   # lower bound on CM about CG
CM_MAX =  0.08   # upper bound on CM about CG
AR_MIN = 4.0
AR_MAX = 12.0
TIP_CHORD_MIN = 0.05

LOG_CSV = "optimization_results.csv"
LOG_FIELDS = [
    "run_id",
    "root_chord", "taper", "sweep", "twist", "span",
    "LD", "CL", "CD", "CM_cg",
    "AR", "wetted_area", "lift", "stall_speed", "rbm", "x_cg",
    "score", "ld_term", "wet_term", "vs_term", "rbm_term",
]

def main():   
    # Define the algorithm
    algorithm = DE(
        pop_size=50,
        variant="DE/rand/1/bin",
        CR=0.9,                  
        F=0.8,                    # This acts as the base mutation weight
        dither="vector",
        jitter=False
    )
    termination = DefaultSingleObjectiveTermination(
        xtol=1e-8,
        cvtol=1e-6,
        ftol=1e-6,
        period=20,
        n_max_gen=100,
        n_max_evals=100000
    )
    problem = DeltaWingProblem()

    init_log()
    print("Starting Optimization...")
    res = minimize(problem, algorithm, termination, seed=1, save_history=True, verbose=True)

    print(f"Optimization finished.")

    # Look up values from the log
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
    wetted      = float(row["wetted_area"])
    vs          = float(row["stall_speed"])
    rbm         = float(row["rbm"])
    cm          = float(row["CM_cg"])
    ar          = float(row["AR"])
    lift        = float(row["lift"])
    x_cg        = float(row["x_cg"])
    score       = float(row["score"])

    _, breakdown = compute_score(ld, wetted, vs, rbm)

    print("\n--- OPTIMAL WING GEOMETRY ---")
    print(f"Score      : {score:.4f}  (weighted sum — higher is better)")
    print(f"\n--- SCORE BREAKDOWN ---")
    for term, contribution in breakdown.items():
        w, ref = SCORING[term]
        print(f"  {term:<16}: contribution={contribution:+.4f}  (weight={w}, reference={ref})")
    print(f"\n--- AERODYNAMICS ---")
    print(f"L/D         : {ld:.4f}")
    print(f"CL          : {float(row['CL']):.4f}")
    print(f"CD          : {float(row['CD']):.4f}")
    print(f"CM_cg       : {cm:.4f}")
    print(f"Total lift  : {lift:.4f} N")
    print(f"\n--- GEOMETRY ---")
    print(f"Wetted area : {wetted:.4f} m²")
    print(f"Aspect ratio: {ar:.4f}  (min: {AR_MIN}  max: {AR_MAX})")
    print(f"Stall speed : {vs:.4f} m/s")
    print(f"Root BM     : {rbm:.4f} N·m")
    print(f"CG          : {x_cg:.4f} m  (aft of root LE,  SM={STATIC_MARGIN*100:.0f}% MAC)")
    print(f"\n--- PARAMETERS ---")
    print(f"Params      : Root={best_root:.4f}  Taper={best_taper:.4f}  Sweep={best_sweep:.4f}  Twist={best_twist:.4f}  Span={best_span:.4f}")
    print(f"\nFeasible designs logged to: {LOG_CSV}")
    print(" ")
    
    stl_path, _ = generate_wing("Optimized_Wing", best_span, best_root, best_taper, best_sweep, 0.0, best_twist, airfoil_file)
    visualize_stl(stl_path)
        
class DeltaWingProblem(ElementwiseProblem):
    def __init__(self):
        super().__init__(
            n_var=5,             # Number of variables
            n_obj=1,             # Number of objectives
            n_constr=4,          # Number of constraints
            xl=np.array([0.1, 0.05, 0.0, -15.0, 0.5]), # Lower bounds for variables
            xu=np.array([0.25, 1.0, 60.0, 5.0, 0.9])   # Upper bounds for variables
        )

    def _evaluate(self, x, out, *args, **kwargs):
        root_chord = x[0]
        taper      = x[1]
        sweep      = x[2]
        twist      = x[3]
        span       = x[4]

        run_id = f"wing_{uuid.uuid4().hex[:8]}" 
        
        try:
            stl_path, analysis_path = generate_wing(run_id, span, root_chord, taper, sweep, 0.0, twist, airfoil_file)
            Sref = 0.5 * (root_chord + root_chord * taper) * span
            x_cg = calc_cg(root_chord, taper, span, sweep)
            aero  = vsp_point(analysis_path, velocity, alpha, Sref, span, root_chord, x_cg)
            CL    = aero["CL"]
            CD    = aero["CD"]
            LD    = aero["LD"]
            CM_cg = aero["CM"]
            lift  = 0.5 * CL * 1.225 * velocity * velocity * Sref

            AR          = aspect_ratio(root_chord, taper, span)
            wetted_area = estimate_wetted_area(root_chord, taper, span)
            vs          = stall_speed(CL, Sref)
            rbm         = root_bending_moment(root_chord, taper, span, CL)
            score, breakdown = compute_score(LD, wetted_area, vs, rbm)

            print(
                f"  [{run_id}] L/D={LD:.3f}  wetted={wetted_area:.4f} m²  "
                f"CM_cg={CM_cg:.4f}  AR={AR:.2f}  Vs={vs:.2f} m/s  RBM={rbm:.3f} N·m  "
                + "  ".join(f"{k}={v:+.4f}" for k, v in breakdown.items())
                + f"  → score={score:.4f}"
            )

            out["F"] = [-score]
            g_lift = MAX_WEIGHT - lift
            g_cm   = max(CM_cg - CM_MAX, CM_MIN - CM_cg)
            g_ar   = max(AR_MIN - AR, AR - AR_MAX)
            g_tip  = TIP_CHORD_MIN - (root_chord * taper)
            out["G"] = [g_lift, g_cm, g_ar, g_tip]

            # Only log feasible designs (all constraints satisfied)
            if g_lift <= 0 and g_cm <= 0 and g_ar <= 0 and g_tip <= 0:
                x_cg_val = calc_cg(root_chord, taper, span, sweep)
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
                    "wetted_area" : round(wetted_area,        6),
                    "lift"        : round(lift,               6),
                    "stall_speed" : round(vs,                 6),
                    "rbm"         : round(rbm,                6),
                    "x_cg"        : round(x_cg_val,          6),
                    "score"       : round(score,              6),
                    "ld_term"     : round(breakdown["ld_ratio"],       6),
                    "wet_term"    : round(breakdown["wetted_area"],     6),
                    "vs_term"     : round(breakdown["stall_speed"],     6),
                    "rbm_term"    : round(breakdown["bending_moment"],  6),
                })

        except Exception as e:
            print(f"Run {run_id} failed: {e}")
            out["F"] = [1e10]
            out["G"] = [1e10, 1e10, 1e10, 1e10]
        finally:
            for filename in glob.glob(f"{run_id}*"):
                try:
                    os.remove(filename)
                except OSError:
                    pass

def generate_wing(wing_name, wingspan, root_chord, taper_ratio, sweep_angle, dihedral_angle, twist_angle, airfoil_file):
    tip_chord = root_chord * taper_ratio
    airfoil_fwd = airfoil_file.replace("\\", "/")

    script_lines = [
        "void main() {",
        "    VSPCheckSetup();",
        "    ClearVSPModel();",
        f'    string wing_id = AddGeom( "WING" );',
        f'    SetParmVal( wing_id, "TotalSpan",      "WingGeom", {wingspan} );',
        f'    SetParmVal( wing_id, "Root_Chord",     "XSec_1",   {root_chord} );',
        f'    SetParmVal( wing_id, "Tip_Chord",      "XSec_1",   {tip_chord} );',
        f'    SetParmVal( wing_id, "Sweep",          "XSec_1",   {sweep_angle} );',
        f'    SetParmVal( wing_id, "Dihedral",       "XSec_1",   {dihedral_angle} );',
        f'    SetParmVal( wing_id, "Twist",          "XSec_1",   {twist_angle} );',
        f'    SetParmVal( wing_id, "Twist_Location", "XSec_1",   1.0 );',
        f'    SetParmVal( wing_id, "SectTess_U",     "XSec_1",   {wing_span_res}.0 );',
        f'    SetParmVal( wing_id, "Tess_W",         "Shape",    {wing_chord_res}.0 );',
        f'    string root_surf = GetXSecSurf( wing_id, 0 );',
        f'    ChangeXSecShape( root_surf, 0, XS_FILE_AIRFOIL );',
        f'    string root_xsec = GetXSec( root_surf, 0 );',
        f'    ReadFileAirfoil( root_xsec, "{airfoil_fwd}" );',
        f'    string tip_surf = GetXSecSurf( wing_id, 1 );',
        f'    ChangeXSecShape( tip_surf, 1, XS_FILE_AIRFOIL );',
        f'    string tip_xsec = GetXSec( tip_surf, 1 );',
        f'    ReadFileAirfoil( tip_xsec, "{airfoil_fwd}" );',
        f'    SetSetFlag( wing_id, 1, true );',
        f'    Update();',
        f'    WriteVSPFile( "{wing_name}.vsp3", SET_ALL );',
        f'    ExportFile( "{wing_name}.stl", 0, EXPORT_STL );',
        "}",
    ]

    script_path = f"{wing_name}_geom.vspscript"
    with open(script_path, 'w') as f:
        f.write("\n".join(script_lines))

    print(f"--- Running geometry generation ({script_path}) ---")
    subprocess.run([vsp_exe, "-script", script_path], check=True)
    os.remove(script_path)

    stl_path = f"{wing_name}.stl"
    vsp3_path = f"{wing_name}.vsp3"
    print(f"STL generated: {stl_path}")
    print(f"VSP file saved: {vsp3_path}")
    return stl_path, vsp3_path

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

    script_lines = [
        "void main() {",
        f'    ClearVSPModel();',
        f'    ReadVSPFile( "{vsp3_path}" );',
        f'    SetAnalysisInputDefaults( "VSPAEROComputeGeometry" );',
        f'    array< int > thick_set = GetIntAnalysisInput( "VSPAEROComputeGeometry", "GeomSet" );',
        f'    array< int > thin_set = GetIntAnalysisInput( "VSPAEROComputeGeometry", "ThinGeomSet" );',
        f'    thick_set[0] = ( SET_TYPE::SET_NONE );',
        f'    thin_set[0] = ( SET_TYPE::SET_ALL );',
        f'    SetIntAnalysisInput( "VSPAEROComputeGeometry", "GeomSet", thick_set );',
        f'    SetIntAnalysisInput( "VSPAEROComputeGeometry", "ThinGeomSet", thin_set );',
        f'    Print( "--- Running Meshing ---" );',
        f'    ExecAnalysis( "VSPAEROComputeGeometry" );',
        f'    SetAnalysisInputDefaults( "VSPAEROSweep" );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Sref",           {darr(Sref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "cref",           {darr(cref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "bref",           {darr(bref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Xcg",            {darr(x_cg)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaStart",     {darr(float(alpha))}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaEnd",       {darr(float(alpha))}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "AlphaNpts",      {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "MachStart",      {darr(mach)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "MachNpts",       {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Vinf",           {darr(100.0)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "WakeNumIter",    {iarr(15)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "NCPU",           {iarr(8)}, 0 );',
        f'    Print( "--- Running Aero Point ---" );',
        f'    ExecAnalysis( "VSPAEROSweep" );',
        "}",
    ]

    script_path = f"{vsp3_path.replace('.vsp3', '')}_aero.vspscript"
    with open(script_path, 'w') as f:
        f.write("\n".join(script_lines))

    print(f"--- Running aero point (alpha={alpha}) ---")
    subprocess.run([vsp_exe, "-script", script_path], check=True)
    os.remove(script_path)

    polar_file = vsp3_path.replace(".vsp3", ".polar")
    CL, CD, Cm = parse_polar(polar_file)
    cl = CL[0]
    cd = CD[0]
    cm = Cm[0]
    return {
        "CL": cl,
        "CD": cd,
        "LD": cl / cd,
        "CM": cm,
    }

def parse_polar(polar_path):
    CL, CD, Cm = [], [], []
    col_cl = col_cd = col_cm = None
 
    with open(polar_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            tokens = stripped.split()
            if tokens[0] == 'Beta':
                col_cl = tokens.index('CLtot')
                col_cd = tokens.index('CDtot')
                col_cm = tokens.index('CMytot')
                continue
            if col_cl is None:
                continue
            try:
                CL.append(float(tokens[col_cl]))
                CD.append(float(tokens[col_cd]))
                Cm.append(float(tokens[col_cm]))
            except (ValueError, IndexError):
                continue
 
    return CL, CD, Cm

def compute_score(ld, wetted_area, stall_spd, rbm):
    w_ld,  ref_ld  = SCORING["ld_ratio"]
    w_wet, ref_wet = SCORING["wetted_area"]
    w_vs,  ref_vs  = SCORING["stall_speed"]
    w_rbm, ref_rbm = SCORING["bending_moment"]

    term_ld  =  w_ld  * (ld         / ref_ld)   # positive: reward high L/D
    term_wet = -w_wet * (wetted_area / ref_wet)  # negative: penalise large area
    term_vs  = -w_vs  * (stall_spd  / ref_vs)   # negative: penalise high stall speed
    term_rbm = -w_rbm * (rbm        / ref_rbm)  # negative: penalise high bending moment

    score = term_ld + term_wet + term_vs + term_rbm

    breakdown = {
        "ld_ratio"      : term_ld,
        "wetted_area"   : term_wet,
        "stall_speed"   : term_vs,
        "bending_moment": term_rbm,
    }
    return score, breakdown

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
            if best is None or float(row["score"]) > float(best["score"]):
                best = row
    return best

# AngelScript array helpers
def iarr(v):
    return f"array<int> = {{{v}}}"

def darr(v):
    return f"array<double> = {{{v}}}"

if __name__ == "__main__":
    main()