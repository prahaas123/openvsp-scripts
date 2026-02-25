import os
import openvsp as vsp
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

wing_span_res = 20
wing_chord_res = 50
velocity = 10 # m/s
alpha = 5 # degrees AoA

airfoil_file = r"C:\Users\kprah\Desktop\Prahaas\WatArrow\CFD Automation\Airfoils\dae21.dat"

def main():    
    # Define the algorithm
    algorithm = NSGA2(pop_size=200, eliminate_duplicates=True)
    termination = get_termination("n_gen", 50)
    problem = DeltaWingProblem()

    print("Starting Optimization...")
    res = minimize(problem, algorithm, termination, seed=1, save_history=True, verbose=True)

    print(f"Optimization finished. Found {len(res.X)} solutions in the Pareto front.")
    
    final_objectives = -res.F 
    
    data = np.hstack([res.X, final_objectives])
    columns = ["Root_Chord", "Taper", "Dihedral", "Twist", "Span", "L_D", "Lift"] 
    df = pd.DataFrame(data, columns=columns)
    df.to_csv("optimization_results.csv", index=False)
    print("--- OPTIMIZATION COMPLETE. Saved to optimization_results.csv ---")

    for i in range(len(res.X)):
        print(f"Design {i}: L/D={final_objectives[i,0]:.4f}, Lift={final_objectives[i,1]:.4f}")
        print(f"   Params: Root={res.X[i,0]:.2f}, Taper={res.X[i,1]:.2f}, Dihedral={res.X[i,2]:.2f}, Twist={res.X[i,3]:.2f}, Span={res.X[i,4]:.2f}")
        
    plot_parallel_coordinates()
        
class DeltaWingProblem(ElementwiseProblem):
    def __init__(self):
        super().__init__(
            n_var=5,             # Number of variables 
            n_obj=2,             # Number of objectives 
            n_constr=0,          # Constraints 
            xl=np.array([50, 0.05, -10.0, -15.0, 50]), # Lower bounds for variables
            xu=np.array([200, 1.0, 10.0, 15.0, 200])  # Upper bounds for variables
        )

    def _evaluate(self, x, out, *args, **kwargs):
        root_chord = x[0]
        taper      = x[1]
        dihedral   = x[2]
        twist      = x[3]
        span       = x[4]
        
        airfoil = r"C:\Users\kprah\Desktop\Prahaas\WatArrow\CFD Automation\Airfoils\dae21.dat"

        run_id = f"wing_{uuid.uuid4().hex[:8]}" 
        
        try:
            stl_path, analysis_path = generate_wing(run_id, span, root_chord, taper, 0.0, dihedral, twist, airfoil)
            CL, CD, LD = vsp_point(analysis_path, velocity, alpha, 0.5 * (root_chord + root_chord * taper) * span, span, root_chord)
            lift = 2 * CL * 1.225 * velocity * velocity * root_chord * span / 10000
            out["F"] = [-LD, -lift]
        except Exception as e:
            print(f"Run {run_id} failed: {e}")
            out["F"] = [1e10, 1e10] # Huge penalty
        finally:
            for filename in glob.glob(f"{run_id}*"):
                try:
                    os.remove(filename)
                except OSError:
                    pass

def generate_wing(wing_name, wingspan, root_chord, taper_ratio, sweep_angle, dihedral_angle, twist_angle, airfoil_file):
    vsp.VSPCheckSetup()
    vsp.ClearVSPModel()

    # Create the wing
    wing_id = vsp.AddGeom("WING")
    
    # Setting parameters
    tip_chord = root_chord * taper_ratio
    
    vsp.SetParmVal(wing_id, "TotalSpan", "WingGeom", wingspan)
    vsp.SetParmVal(wing_id, "Root_Chord", "XSec_1", root_chord)
    vsp.SetParmVal(wing_id, "Tip_Chord", "XSec_1", tip_chord)
    vsp.SetParmVal(wing_id, "Sweep", "XSec_1", sweep_angle)
    vsp.SetParmVal(wing_id, "Dihedral", "XSec_1", dihedral_angle)
    vsp.SetParmVal(wing_id, "Twist", "XSec_1", twist_angle)
    vsp.SetParmVal(wing_id, "Twist_Location", "XSec_1", 1.0)
    vsp.SetParmVal(wing_id, "SectTess_U", "XSec_1", wing_span_res)
    vsp.SetParmVal(wing_id, "Tess_W", "Shape", wing_chord_res)
    
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

    # Finalize and xxport
    vsp.Update()
    
    stl_path = f"{wing_name}.stl"
    analysis_path = f"{wing_name}.vsp3"
    vsp.WriteVSPFile(analysis_path)
    vsp.ExportFile(stl_path, 0, vsp.EXPORT_STL)
    print(f"STL generated: {stl_path}")
    print(f"VSP file saved: {analysis_path}")
    
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

def vsp_point(fname_vspaerotests, vin, alpha, Sref, bref, cref):
    # Load model
    vsp.ReadVSPFile(fname_vspaerotests)
    
    # Geometry
    geom_analysis = "VSPAEROComputeGeometry"
    vsp.SetAnalysisInputDefaults(geom_analysis)
    
    vsp.SetIntAnalysisInput(geom_analysis, "GeomSet", [vsp.SET_NONE])      
    vsp.SetIntAnalysisInput(geom_analysis, "ThinGeomSet", [vsp.SET_SHOWN])
    
    print(f"--- Running Meshing ({geom_analysis}) ---")
    vsp.ExecAnalysis(geom_analysis)

    aero_analysis = "VSPAEROSweep"
    vsp.SetAnalysisInputDefaults(aero_analysis)
    
    # Reference dimensions
    vsp.SetDoubleAnalysisInput(aero_analysis, "Sref", [Sref])
    vsp.SetDoubleAnalysisInput(aero_analysis, "cref", [cref])
    vsp.SetDoubleAnalysisInput(aero_analysis, "bref", [bref])
    
    # Flight conditions    
    vsp.SetDoubleAnalysisInput(aero_analysis, "AlphaStart", [alpha])
    vsp.SetIntAnalysisInput(aero_analysis, "AlphaNpts", [1])    
    mach = vin / 343.0
    vsp.SetDoubleAnalysisInput(aero_analysis, "MachStart", [mach])
    vsp.SetIntAnalysisInput(aero_analysis, "MachNpts", [1])
    vsp.SetIntAnalysisInput(aero_analysis, "WakeNumIter", [15]) 
    vsp.SetDoubleAnalysisInput(aero_analysis, "Vinf", [100.0])

    print(f"--- Running Aero Sweep ({aero_analysis}) ---")
    rid = vsp.ExecAnalysis(aero_analysis)

    # Results
    polar_res = vsp.FindLatestResultsID("VSPAERO_Polar")
    cl = vsp.GetDoubleResults(polar_res, "CLtot")[0]
    cd = vsp.GetDoubleResults(polar_res, "CDtot")[0]
    return cl, cd, cl/cd

def plot_parallel_coordinates():
    df = pd.read_csv("optimization_results.csv")
    fig = px.parallel_coordinates(
        df,
        color="Lift",
        color_continuous_scale='Plasma',
        dimensions=["Root_Chord", "Taper Ratio", "Sweep", "Span"],
        title="Design Genetic Optimization: Geometry vs Lift",
        template="plotly_dark"
    )
    fig.show()
    fig2= px.parallel_coordinates(
        df,
        color="L_D",
        color_continuous_scale='Plasma',
        dimensions=["Root_Chord", "Taper Ratio", "Sweep", "Span"],
        title="Design Genetic Optimization: Geometry vs L/D",
        template="plotly_dark"
    )
    fig2.show()

if __name__ == "__main__":
    plot_parallel_coordinates()