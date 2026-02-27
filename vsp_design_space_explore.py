import os
import openvsp as vsp
import pyvista as pv
import numpy as np
import uuid
import glob
import pandas as pd
import itertools
import os
import time
from tqdm import tqdm
import plotly.express as px
import plotly.graph_objects as go

wing_span_res = 20
wing_chord_res = 50
velocity = 10 # m/s
alpha = 5 # degrees AoA

airfoil_file = r"C:\Users\kprah\Desktop\Prahaas\WatArrow\CFD Automation\Airfoils\dae21.dat"

ranges = {
    "root_chord": np.linspace(0.5, 2.0, 4),     # meters
    "taper":      np.linspace(0.1, 0.4, 4),     # ratio
    "sweep":      np.linspace(40.0, 70.0, 4),   # degrees
    "twist":      np.linspace(-10.0, 0.0, 3),   # degrees
    "span":       np.linspace(0.5, 2.0, 3)      # meters
}

def main():
    keys = list(ranges.keys())
    values = list(ranges.values())
    grid_points = list(itertools.product(*values))
    
    total_sims = len(grid_points)
    print(f"--- STARTING PARAMETRIC SWEEP ---")
    print(f"Total Designs to Simulate: {total_sims}")
    print(f"Estimated Time (@5s/sim): {total_sims * 5 / 60:.1f} minutes")
    time.sleep(5)  # Brief pause before starting
    
    results = []

    for point in tqdm(grid_points):
        r_c, tap, swe, twi, wing = point
        run_id = f"wing_{uuid.uuid4().hex[:8]}"
        
        try:
            stl_path, analysis_path = generate_wing(run_id, wing, r_c, tap, swe, 0.0, twi, airfoil_file)
            CL, CD, LD = vsp_point(analysis_path, velocity, alpha, 0.5 * (r_c + r_c * tap) * wing, wing, r_c)
            
            results.append({
                "Root Chord": r_c,
                "Taper": tap,
                "Sweep": swe,
                "Twist": twi,
                "Wingspan": wing,
                "CL": CL,
                "CD": CD,
                "L_D": LD,
                "Lift": 2 * CL * 1.225 * velocity * velocity * r_c * wing
            })
            
        except Exception as e:
            print(f"Failed on design {point}: {e}")
            
        for filename in glob.glob(f"{run_id}*"):
            try:
                os.remove(filename)
            except OSError:
                pass

    df = pd.DataFrame(results)
    df.to_csv("sweep_results.csv", index=False)
    print("--- SWEEP COMPLETE. Saved to sweep_results.csv ---")
    
    return df

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
    vsp.SetIntAnalysisInput(aero_analysis, "NCPU", [8])

    print(f"--- Running Aero Sweep ({aero_analysis}) ---")
    rid = vsp.ExecAnalysis(aero_analysis)

    # Results
    polar_res = vsp.FindLatestResultsID("VSPAERO_Polar")
    cl = vsp.GetDoubleResults(polar_res, "CLtot")[0]
    cd = vsp.GetDoubleResults(polar_res, "CDtot")[0]
    return cl, cd, cl/cd

def plot_parallel_coordinates():
    df = pd.read_csv("sweep_results.csv")
    threshold = df["L_D"].quantile(0.95) # Top 5% L/D as threshold
    df_top = df[df["L_D"] >= threshold]
    fig = px.parallel_coordinates(
        df_top,
        dimensions=["Root Chord", "Taper", "Sweep", "Twist", "Wingspan"],
        color="L_D",
        color_continuous_scale="Plasma",
        title=f"Design Space: Top 5% of L/D Results (L/D > {threshold:.2f})",
        template="plotly_dark"
    )
    fig.show()
    
def plot_pareto_front():
    df = pd.read_csv("sweep_results.csv")
    df_sorted = df.sort_values(by='CL', ascending=False)
    pareto_front = []
    current_max_ld = -1.0
    
    for _, row in df_sorted.iterrows():
        if row['L_D'] > current_max_ld:
            pareto_front.append(row)
            current_max_ld = row['L_D']
            
    pareto_df = pd.DataFrame(pareto_front)
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df['CL'], y=df['L_D'],
        mode='markers',
        name='All Designs',
        marker=dict(color='lightgrey', size=5, opacity=0.3)
    ))
    fig.add_trace(go.Scatter(
        x=pareto_df['CL'], y=pareto_df['L_D'],
        mode='lines+markers',
        name='Pareto Front (Optimal)',
        marker=dict(color='yellow', size=8),
        line=dict(dash='dash')
    ))
    fig.update_layout(
        title="Pareto Front: CL vs L/D Optimization",
        xaxis_title="Lift Coefficient (CL)",
        yaxis_title="Lift-to-Drag Ratio (L/D)",
        template="plotly_dark"
    )
    fig.show()
    
def plot_splom():
    df = pd.read_csv("sweep_results.csv")
    fig = px.scatter_matrix(
        df,
        dimensions=["Root Chord", "Taper", "Sweep", "Twist", "Wingspan"],
        color="L_D",
        color_continuous_scale="Viridis",
        title="Pairwise Relationships in Design Space",
        template="plotly_dark"
    )
    fig.update_traces(diagonal_visible=False) # Hides the diagonal plots
    fig.update_layout(height=800, width=1000)
    fig.show()
    
if __name__ == "__main__":
    # main()
    plot_parallel_coordinates()
    plot_pareto_front()
    plot_splom()