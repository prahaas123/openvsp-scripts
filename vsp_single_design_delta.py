import os
import pyvista as pv
import numpy as np
import pandas as pd
import csv
import glob
import subprocess
import plotly.graph_objects as go
import plotly.express as px
from dash import Dash, dcc, html, Input, Output
import dash_bootstrap_components as dbc

vsp_exe = r"C:\Program Files\OpenVSP-3.47.0\vsp.exe"

wing_span_res = 20
wing_chord_res = 50
velocities = list(range(10.0, 50.0, 5.0)) # m/s
alphas = list(range(-5, 15, 1)) # degrees AoA

airfoil_file = r"Airfoils\mh45.dat"

root_chord = 0.232
taper_ratio = 0.4288
sweep = 27.2141
dihedral = 3.0
twist = -0.1073
span = 1.0589
x_cg = 0.1529
elevon_length = 0.35        # % of chord
elevon_start = 0.4          # % of wingspan
elevon_end = 0.9            # % of wingspan

def main():
    try:
        os.remove("aero_full.csv")
        os.remove("stability.csv")
    except OSError:
        pass
    stl_path, vsp3_path = generate_wing("wing")
    visualize_stl(stl_path)
    Sref = 0.5 * (root_chord + root_chord * taper_ratio) * span
    bref = span
    cref = root_chord

    # Aero sweep
    csv_exists = False
    for v in velocities:
        print(f"\n=== Running VSP Aero Sweep at {v} m/s ===")
        CL, CD, CDi, Cm = vsp_sweep(vsp3_path, v, Sref, bref, cref)
        lod_filename = "wing.lod"
        cl_data = read_lift_distribution(lod_filename)
        first_aoa_key = list(cl_data.keys())[0]
        span_locations = cl_data[first_aoa_key]['SoverB']
        aero_headers = ["Velocity", "Alpha_deg", "CL", "CD", "Cm"] + [f"Cl_span_{loc:.4f}" for loc in span_locations] + ["Oswald_efficiency"]

        # Combine sweep results with local spanwise results
        aero_results = []
        for i, alpha in enumerate(alphas):
            e = [compute_oswald(CL[i], CDi[i], Sref, bref)]
            base_row = [v, alpha, CL[i], CD[i], Cm[i]]
            # Safely find the corresponding AoA in the parsed dictionary to avoid errors
            closest_aoa_key = min(cl_data.keys(), key=lambda k: abs(k - alpha))
            spanwise_cls = cl_data[closest_aoa_key]['Cl']
            full_row = base_row + spanwise_cls + e
            aero_results.append(full_row)

        # Write to CSV
        aero_filename = "aero_full.csv"
        with open(aero_filename, 'a', newline='') as f:
            writer = csv.writer(f)
            if not csv_exists:
                writer.writerow(aero_headers)
                csv_exists = True
            writer.writerows(aero_results)

        # File Cleanup
        for filename in glob.glob("wing*"):
            try:
                os.remove(filename)
            except OSError:
                pass

    # Stability sweep
    stab_csv_exists = False
    
    for v in velocities:
        stl_path, vsp3_path = generate_wing("wing")
        vsp_stability(vsp3_path, v, Sref, bref, cref)
        stab_dict = read_stability("wing.stab")
        sm = get_sm("wing.stab", cg=x_cg, mac=cref)
        
        stability_filename = "stability.csv"
        stab_headers = ["Velocity"] + list(stab_dict.keys()) + ["StaticMargin"]
        stab_columns = [v] + list(stab_dict.values()) + [sm]
        
        with open(stability_filename, 'a', newline='') as f:
            writer = csv.writer(f)
            if not stab_csv_exists:
                writer.writerow(stab_headers)
                stab_csv_exists = True
            writer.writerow(stab_columns)

        for filename in glob.glob("wing*"):
            try:
                os.remove(filename)
            except OSError:
                pass

def generate_wing(wing_name):
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
        f'    AddSubSurf( wing_id, SS_CONTROL );',
        f'    SetParmVal( wing_id, "SE_Const_Flag",  "SS_Control_1", 1.0 );',
        f'    SetParmVal( wing_id, "Length_C_Start", "SS_Control_1", {elevon_length} );',
        f'    SetParmVal( wing_id, "EtaFlag",        "SS_Control_1", 1.0 );',
        f'    SetParmVal( wing_id, "EtaStart",       "SS_Control_1", {elevon_start} );',
        f'    SetParmVal( wing_id, "EtaEnd",         "SS_Control_1", {elevon_end} );',
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

def vsp_sweep(vsp3_path, v, Sref, bref, cref):
    mach = v / 343.0

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
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Sref",        {darr(Sref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "cref",        {darr(cref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "bref",        {darr(bref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaStart",  {darr(float(alphas[0]))}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaEnd",    {darr(float(alphas[-1]))}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "AlphaNpts",   {iarr(len(alphas))}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "MachStart",   {darr(float(mach))}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "MachNpts",    {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Vinf",        {darr(v)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Rho",         {darr(1.225)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Xcg",         {darr(x_cg)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "WakeNumIter", {iarr(15)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "NCPU",        {iarr(8)}, 0 );',
        f'    Print( "--- Running Aero Sweep ---" );',
        f'    SetStringAnalysisInput( "VSPAEROSweep", "RedirectFile", array<string> = {{"{vsp3_path}_log.txt"}}, 0 );',
        f'    ExecAnalysis( "VSPAEROSweep" );',
        "}",
    ]

    script_path = "wing_sweep.vspscript"
    with open(script_path, 'w') as f:
        f.write("\n".join(script_lines))

    print(f"--- Running aero sweep ({script_path}) ---")
    subprocess.run([vsp_exe, "-script", script_path], check=True)
    os.remove(script_path)

    CL, CD, CDi, Cm = parse_polar("wing.polar")
    return CL, CD, CDi, Cm

def vsp_stability(vsp3_path, v, Sref, bref, cref):
    mach = v / 343.0

    script_lines = [
        "void main() {",
        f'    ReadVSPFile( "{vsp3_path}" );',
        f'    SetAnalysisInputDefaults( "VSPAEROComputeGeometry" );',
        f'    array< int > thick_set = GetIntAnalysisInput( "VSPAEROComputeGeometry", "GeomSet" );',
        f'    array< int > thin_set = GetIntAnalysisInput( "VSPAEROComputeGeometry", "ThinGeomSet" );',
        f'    thick_set[0] = ( SET_TYPE::SET_NONE );',
        f'    thin_set[0] = ( SET_TYPE::SET_ALL );',
        f'    int cs_pitch_id = CreateVSPAEROControlSurfaceGroup();',
        f'    SetVSPAEROControlGroupName( "Pitch", cs_pitch_id );',
        f'    AddAllToVSPAEROControlSurfaceGroup( cs_pitch_id );',
        f'    int cs_roll_id = CreateVSPAEROControlSurfaceGroup();',
        f'    SetVSPAEROControlGroupName( "Roll", cs_roll_id );',
        f'    AddAllToVSPAEROControlSurfaceGroup( cs_roll_id );',
        f'    Update();',
        f'    string container_id = FindContainer( "VSPAEROSettings", 0 );',
        f'    array<string> geoms = FindGeoms();',
        f'    string wid = geoms[0];',
        f'    string cs_id = GetSubSurf( wid, 0 );',
        f'    string group_str = "ControlSurfaceGroup_" + formatInt( cs_pitch_id + 1, "" );',
        f'    string parm_id = FindParm( container_id, "Surf_" + cs_id + "_1_Gain", group_str );',
        f'    SetParmVal( parm_id, -1.0 );',
        f'    Print( "--- Running Meshing (stability) ---" );',
        f'    ExecAnalysis( "VSPAEROComputeGeometry" );',
        f'    SetAnalysisInputDefaults( "VSPAEROSweep" );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Sref",         {darr(Sref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "cref",         {darr(cref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "bref",         {darr(bref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaStart",   {darr(0.0)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaEnd",     {darr(0.0)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "AlphaNpts",    {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "MachStart",    {darr(mach)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "MachNpts",     {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Vinf",         {darr(100.0)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Xcg",          {darr(x_cg)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "WakeNumIter",  {iarr(15)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "UnsteadyType", {iarr(1)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "NCPU",         {iarr(8)}, 0 );',
        f'    Print( "--- Running Stability Sweep ---" );',
        f'    SetStringAnalysisInput( "VSPAEROSweep", "RedirectFile", array<string> = {{"{vsp3_path}_log.txt"}}, 0 );',
        f'    ExecAnalysis( "VSPAEROSweep" );',
        "}",
    ]

    script_path = "wing_stab.vspscript"
    with open(script_path, 'w') as f:
        f.write("\n".join(script_lines))

    print(f"--- Running stability sweep ({script_path}) ---")
    subprocess.run([vsp_exe, "-script", script_path], check=True)
    os.remove(script_path)

def parse_polar(polar_path):
    CL, CD, CDi, Cm = [], [], [], []
    col_cl = col_cd = col_cdi = col_cm = None
 
    with open(polar_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            tokens = stripped.split()
            if tokens[0] == 'Beta':
                col_cl   = tokens.index('CLtot')
                col_cd   = tokens.index('CDtot')
                col_cdi  = tokens.index('CDi')
                col_cm   = tokens.index('CMytot')
                continue
            if col_cl is None:
                continue
            try:
                CL.append(float(tokens[col_cl]))
                CD.append(float(tokens[col_cd]))
                CDi.append(float(tokens[col_cdi]))
                Cm.append(float(tokens[col_cm]))
            except (ValueError, IndexError):
                continue
 
    return CL, CD, CDi, Cm

def read_stability(stab_path, output_file="stability.csv"):
    col_alpha = col_beta = col_p = col_q = col_r = col_pitch = col_roll = None
    rows = {}
 
    with open(stab_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            tokens = stripped.split()
            if tokens[0] == 'Coef':
                col_alpha = tokens.index('Alpha')
                col_beta  = tokens.index('Beta')
                col_p     = tokens.index('p')
                col_q     = tokens.index('q')
                col_r     = tokens.index('r')
                col_pitch = tokens.index('ConGrp_1')
                col_roll  = tokens.index('ConGrp_2')
                continue
            if col_alpha is None:
                continue
            coef = tokens[0]
            if coef in ('CL', 'CS', 'CMl', 'CMm', 'CMn', 'CD'):
                try:
                    rows[coef] = [float(tokens[col_alpha]),
                                  float(tokens[col_beta]),
                                  float(tokens[col_p]),
                                  float(tokens[col_q]),
                                  float(tokens[col_r]),
                                  float(tokens[col_pitch]),
                                  float(tokens[col_roll])]
                except (ValueError, IndexError):
                    pass
 
    vsp_dict = {}
    if 'CL'  in rows: vsp_dict['CL_de']   = rows['CL'] [5]  # CL   wrt ConGrp_1 (pitch)
    if 'CS'  in rows: vsp_dict['CY_beta'] = rows['CS'] [1]  # CS   wrt Beta
    if 'CS'  in rows: vsp_dict['CY_p']    = rows['CS'] [2]  # CS   wrt p
    if 'CS'  in rows: vsp_dict['CY_r']    = rows['CS'] [4]  # CS   wrt r
    if 'CMl' in rows: vsp_dict['Cl_beta'] = rows['CMl'][1]  # CMl  wrt Beta
    if 'CMl' in rows: vsp_dict['Cl_p']    = rows['CMl'][2]  # CMl  wrt p
    if 'CMl' in rows: vsp_dict['Cl_r']    = rows['CMl'][4]  # CMl  wrt r
    if 'CMl' in rows: vsp_dict['Cl_da']   = rows['CMl'][6]  # CMl  wrt ConGrp_2 (roll)
    if 'CMm' in rows: vsp_dict['Cm_q']    = rows['CMm'][3]  # CMm  wrt q
    if 'CMm' in rows: vsp_dict['Cm_de']   = rows['CMm'][5]  # CMm  wrt ConGrp_1 (pitch)
    if 'CMn' in rows: vsp_dict['Cn_beta'] = rows['CMn'][1]  # CMn  wrt Beta
    if 'CMn' in rows: vsp_dict['Cn_p']    = rows['CMn'][2]  # CMn  wrt p
    if 'CMn' in rows: vsp_dict['Cn_r']    = rows['CMn'][4]  # CMn  wrt r
    if 'CMn' in rows: vsp_dict['Cn_da']   = rows['CMn'][6]  # CMn  wrt ConGrp_2 (roll)
 
    return vsp_dict
    
def read_lift_distribution(filepath, target_vortex_sheet=1):
    data_by_aoa = {}
    current_aoa = None
    
    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith("AoA_"):
                parts = line.split()
                current_aoa = float(parts[1])
                
                if current_aoa not in data_by_aoa:
                    data_by_aoa[current_aoa] = {'SoverB': [], 'Cl': []}
            
            else:
                parts = line.split()
                if len(parts) > 15:
                    try:
                        vortex_sheet = int(parts[1])
                        if vortex_sheet == target_vortex_sheet and current_aoa is not None:
                            soverb = float(parts[7])
                            cl = float(parts[11])
                            data_by_aoa[current_aoa]['SoverB'].append(soverb)
                            data_by_aoa[current_aoa]['Cl'].append(cl)
                            
                    except ValueError:
                        pass
                        
    return data_by_aoa

def compute_oswald(cl, cdi, s, b):
    AR = (b ** 2) / s
    cl_array = np.array(cl, dtype=float)
    cdi_array = np.array(cdi, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        e = (cl_array ** 2) / (np.pi * AR * cdi_array)
        e = np.where(np.isfinite(e), e, 0.0)
    if e.ndim == 0:
        return float(e)
    return e
    
def get_sm(stab_path, cg=x_cg, mac=1.0):
    with open(stab_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            
            if stripped.startswith('X_np'):
                parts = stripped.split()
                try:
                    np_val = float(parts[1])
                    sm = (np_val - cg) / mac
                    return sm
                    
                except (ValueError, IndexError):
                    return float('nan')
                    
    return float('nan')

# AngelScript array helpers
def iarr(v):
    return f"array<int> = {{{v}}}"

def darr(v):
    return f"array<double> = {{{v}}}"

def plot_dashboards(sweep_csv="aero_full.csv", stab_csv="stability.csv"):
    # --- Data Loading & Pre-processing ---
    df_aero = pd.read_csv(sweep_csv)
    df_stab = pd.read_csv(stab_csv)

    # Calculate L/D Ratio
    df_aero['L_D'] = df_aero['CL'] / df_aero['CD']

    # Extract spanwise columns and their numerical positions
    span_cols = [c for c in df_aero.columns if "Cl_span_" in c]
    span_pos = [float(c.split("_")[-1]) for c in span_cols]

    # Initialize Dash App with a Dark Theme
    app = Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])

    # Default plotly template for dark mode
    PLOT_TEMPLATE = "plotly_dark"

    app.layout = dbc.Container([
        html.H1("Aircraft Performance Dashboard", className="mt-4 mb-4 text-center"),
        
        dcc.Tabs([
            # --- TAB 1: AERODYNAMICS ---
            dcc.Tab(label='Aerodynamics', children=[
                # Row 1: CL, CD, Cm
                dbc.Row([
                    dbc.Col(dcc.Graph(id='graph-cl-alpha'), width=4),
                    dbc.Col(dcc.Graph(id='graph-cd-alpha'), width=4),
                    dbc.Col(dcc.Graph(id='graph-cm-alpha'), width=4),
                ], className="mt-4"),
                
                # Row 2: Polar, L/D, Oswald
                dbc.Row([
                    dbc.Col(dcc.Graph(id='graph-cl-cd'), width=4),
                    dbc.Col(dcc.Graph(id='graph-ld-alpha'), width=4),
                    dbc.Col(dcc.Graph(id='graph-oswald-alpha'), width=4),
                ], className="mt-4"),

                # Row 3: Lift Distribution (Full Width)
                dbc.Row([
                    dbc.Col([
                        html.Div([
                            html.Label("Select Angle of Attack (Alpha):", className="me-2"),
                            dcc.Dropdown(
                                id='alpha-dropdown',
                                options=[{'label': f"{a}°", 'value': a} for a in df_aero['Alpha_deg'].unique()],
                                value=df_aero['Alpha_deg'].min(),
                                clearable=False,
                                style={'color': 'black', 'width': '100px', 'display': 'inline-block', 'verticalAlign': 'middle'}
                            )
                        ], className="text-center mt-4"),
                        dcc.Graph(id='graph-span-lift'),
                    ], width=12)
                ])
            ], className="p-3"),

            # --- TAB 2: STABILITY ---
            dcc.Tab(label='Stability', children=[
                dbc.Row([
                    # Left: Derivative Bars
                    dbc.Col([
                        html.H4("Static Stability Derivatives", className="mt-3 text-center"),
                        html.Label("Select Velocity (m/s):"),
                        dcc.Slider(
                            id='vel-slider',
                            min=df_stab['Velocity'].min(),
                            max=df_stab['Velocity'].max(),
                            step=None,
                            marks={int(v): {'label': str(int(v)), 'style': {'color': 'white'}} for v in df_stab['Velocity'].unique()},
                            value=df_stab['Velocity'].min()
                        ),
                        dcc.Graph(id='graph-stab-bars'),
                        html.Div(id='static-margin-text', className="lead text-center font-weight-bold mt-2")
                    ], width=6),

                    # Right: Trends
                    dbc.Col([
                        html.H4("Stability Trends", className="mt-3 text-center"),
                        dcc.Graph(id='graph-lon-trends'),
                        dcc.Graph(id='graph-lat-trends')
                    ], width=6)
                ])
            ], className="p-3")
        ])
    ], fluid=True)

    # --- Callbacks ---

    @app.callback(
        [Output('graph-cl-alpha', 'figure'),
         Output('graph-cd-alpha', 'figure'),
         Output('graph-cm-alpha', 'figure'),
         Output('graph-cl-cd', 'figure'),
         Output('graph-ld-alpha', 'figure'),
         Output('graph-oswald-alpha', 'figure')],
        [Input('alpha-dropdown', 'value')]
    )
    def update_aero_plots(_):
        # Helper to create styled line plots
        def create_fig(y_col, title, x_col='Alpha_deg'):
            fig = px.line(df_aero, x=x_col, y=y_col, title=title, markers=True, template=PLOT_TEMPLATE)
            fig.update_layout(margin=dict(l=20, r=20, t=40, b=20))
            return fig

        return (
            create_fig('CL', 'CL vs Alpha'),
            create_fig('CD', 'CD vs Alpha'),
            create_fig('Cm', 'Cm vs Alpha'),
            create_fig('CL', 'Drag Polar (CL vs CD)', x_col='CD'),
            create_fig('L_D', 'L/D vs Alpha'),
            create_fig('Oswald_efficiency', 'Oswald Efficiency vs Alpha')
        )

    @app.callback(
        Output('graph-span-lift', 'figure'),
        [Input('alpha-dropdown', 'value')]
    )
    def update_span_plot(selected_alpha):
        row = df_aero[df_aero['Alpha_deg'] == selected_alpha].iloc[0]
        cl_values = [row[col] for col in span_cols]
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=span_pos, y=cl_values, mode='lines+markers', fill='tozeroy', name='Local Cl'))
        fig.update_layout(
            template=PLOT_TEMPLATE,
            title=f"Spanwise Lift Distribution (Alpha = {selected_alpha}°)",
            xaxis_title="Normalized Span Position (y/b)",
            yaxis_title="Local Cl",
            height=450
        )
        return fig

    @app.callback(
        [Output('graph-stab-bars', 'figure'),
         Output('graph-lon-trends', 'figure'),
         Output('graph-lat-trends', 'figure'),
         Output('static-margin-text', 'children')],
        [Input('vel-slider', 'value')]
    )
    def update_stability(selected_vel):
        row = df_stab[df_stab['Velocity'] == selected_vel].iloc[0]
        
        # Derivatives Bar Chart
        derivs = ['CY_beta', 'Cl_beta', 'Cn_beta', 'Cm_q', 'Cm_de']
        fig_bar = px.bar(x=derivs, y=[row[d] for d in derivs], title="Derivatives", template=PLOT_TEMPLATE)
        
        # Longitudinal Trends
        fig_lon = px.line(df_stab, x='Velocity', y=['Cm_q', 'Cm_de'], title="Longitudinal Trends", template=PLOT_TEMPLATE)
        
        # Lateral Trends
        fig_lat = px.line(df_stab, x='Velocity', y=['Cl_beta', 'Cn_beta', 'Cl_p', 'Cl_r'], title="Lateral Trends", template=PLOT_TEMPLATE)
        
        sm_text = f"Static Margin: {row['StaticMargin']*100:.2f}%"
        
        return fig_bar, fig_lon, fig_lat, sm_text

    return app

if __name__ == '__main__':
    app = plot_dashboards()
    app.run(debug=True)