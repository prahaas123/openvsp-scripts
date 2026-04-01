import os
import pyvista as pv
import numpy as np
import uuid
import glob
import pandas as pd
import os
import subprocess
import time
from tqdm import tqdm
import plotly.express as px
import plotly.graph_objects as go

vsp_exe = r"C:\Program Files\OpenVSP-3.47.0\vsp.exe"

wing_span_res = 20
wing_chord_res = 50
velocity = 10 # m/s
alpha = 5 # degrees AoA

airfoil_file = r"C:\Users\kprah\Desktop\Prahaas\WatArrow\CFD Automation\Airfoils\dae21.dat"

bounds = {
    "wing_area":  (0.13, 0.21),     # meters^2
    "taper":      (0.1, 1.0),       # ratio
    "sweep":      (10.0, 70.0),     # degrees
    "twist":      (-10.0, 0.0),     # degrees
    "span":       (0.5, 2.0)        # meters
}

num_mc_samples = 200

def main():
    total_sims = num_mc_samples
    print(f"--- STARTING MONTE CARLO SWEEP ---")
    print(f"Total Designs to Simulate: {total_sims}")
    print(f"Estimated Time (@5s/sim): {total_sims * 5 / 60:.1f} minutes")
    time.sleep(5)  # Brief pause before starting
    
    results = []

    for i in tqdm(range(total_sims)):
        s  = np.random.uniform(*bounds["wing_area"])
        tap  = np.random.uniform(*bounds["taper"])
        swe  = np.random.uniform(*bounds["sweep"])
        twi  = np.random.uniform(*bounds["twist"])
        wing = np.random.uniform(*bounds["span"])
        r_c = (2 * s) / (wing * (1 + tap))
        
        point = (r_c, tap, swe, twi, wing)
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

def generate_wing(wing_name, span, root_chord, taper_ratio, sweep, dihedral, twist, airfoil_file):
    tip_chord = root_chord * taper_ratio
    airfoil_fwd = airfoil_file.replace("\\", "/")

    script_lines = [
        "void main() {",
        "    VSPCheckSetup();",
        "    ClearVSPModel();",
        f'    string wing_id = AddGeom( "WING" );',
        f'    SetParmVal( wing_id, "TotalSpan",      "WingGeom", {span} );',
        f'    SetParmVal( wing_id, "Root_Chord",     "XSec_1",   {root_chord} );',
        f'    SetParmVal( wing_id, "Tip_Chord",      "XSec_1",   {tip_chord} );',
        f'    SetParmVal( wing_id, "Sweep",          "XSec_1",   {sweep} );',
        f'    SetParmVal( wing_id, "Dihedral",       "XSec_1",   {dihedral} );',
        f'    SetParmVal( wing_id, "Twist",          "XSec_1",   {twist} );',
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
        plotter = pv.Plotter(title="Delta Wing")
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

    script_path = "point.vspscript"
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
    
# AngelScript array helpers
def iarr(v):
    return f"array<int> = {{{v}}}"

def darr(v):
    return f"array<double> = {{{v}}}"
    
if __name__ == "__main__":
    main()
    plot_parallel_coordinates()
    plot_pareto_front()
    plot_splom()