import os
import subprocess
import pyvista as pv
import numpy as np
import glob
import uuid
import pandas as pd
import plotly.express as px
from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.termination import get_termination

vsp_exe = r"C:\Program Files\OpenVSP-3.47.0\vsp.exe"

wing_span_res = 20
wing_chord_res = 50
velocity = 10 # m/s
alpha = 5 # degrees AoA

airfoil_file = r"C:\Users\kprah\Desktop\Prahaas\WatArrow\CFD Automation\Airfoils\dae21.dat"

def main():    
    # Define the algorithm
    algorithm = NSGA2(pop_size=10, eliminate_duplicates=True)
    termination = get_termination("n_gen", 10)
    problem = DeltaWingProblem()

    print("Starting Optimization...")
    res = minimize(problem, algorithm, termination, seed=1, save_history=True, verbose=True)

    print(f"Optimization finished. Found {len(res.X)} solutions in the Pareto front.")
    
    final_objectives = -res.F 
    
    data = np.hstack([res.X, final_objectives])
    columns = ["Root_Chord", "Taper", "Dihedral", "Twist", "Span", "L_D", "CD"] 
    df = pd.DataFrame(data, columns=columns)
    df.to_csv("optimization_results.csv", index=False)
    print("--- OPTIMIZATION COMPLETE. Saved to optimization_results.csv ---")

    for i in range(len(res.X)):
        print(f"Design {i}: L/D={final_objectives[i,0]:.4f}, CD={final_objectives[i,1]:.4f}")
        print(f"   Params: Root={res.X[i,0]:.2f}, Taper={res.X[i,1]:.2f}, Dihedral={res.X[i,2]:.2f}, Twist={res.X[i,3]:.2f}, Span={res.X[i,4]:.2f}")
        
    plot_parallel_coordinates()
        
class DeltaWingProblem(ElementwiseProblem):
    def __init__(self):
        super().__init__(
            n_var=5,             # Number of variables 
            n_obj=2,             # Number of objectives 
            n_constr=1,          # Constraints 
            xl=np.array([0.5, 0.05, -10.0, -15.0, 0.5]), # Lower bounds for variables
            xu=np.array([2.0, 1.0, 10.0, 15.0, 2.0])  # Upper bounds for variables
        )

    def _evaluate(self, x, out, *args, **kwargs):
        root_chord = x[0]
        taper      = x[1]
        dihedral   = x[2]
        twist      = x[3]
        span       = x[4]

        run_id = f"wing_{uuid.uuid4().hex[:8]}" 
        
        try:
            stl_path, analysis_path = generate_wing(run_id, span, root_chord, taper, 0.0, dihedral, twist, airfoil_file)
            CL, CD, LD = vsp_point(analysis_path, velocity, alpha, 0.5 * (root_chord + root_chord * taper) * span, span, root_chord)
            lift = 2 * CL * 1.225 * velocity * velocity * root_chord * span
            out["F"] = [-LD, CD]
            out["G"] = 17 - lift
        except Exception as e:
            print(f"Run {run_id} failed: {e}")
            out["F"] = [1e10, 1e10] # Huge penalty
            out["G"] = 1e10
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

def vsp_point(vsp3_path, vin, alpha, Sref, bref, cref):
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
    CL, CD, _ = parse_polar(polar_file)
    cl = CL[0]
    cd = CD[0]
    return cl, cd, cl / cd

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

def plot_parallel_coordinates():
    df = pd.read_csv("optimization_results.csv")
    fig = px.parallel_coordinates(
        df,
        color="L_D",
        color_continuous_scale='Plasma',
        dimensions=["Root_Chord", "Taper", "Dihedral", "Twist", "Span", "L_D"],
        title="Design Genetic Optimization: Geometry vs L/D",
        template="plotly_dark"
    )
    fig.show()

# AngelScript array helpers
def iarr(v):
    return f"array<int> = {{{v}}}"

def darr(v):
    return f"array<double> = {{{v}}}"

if __name__ == "__main__":
    main()