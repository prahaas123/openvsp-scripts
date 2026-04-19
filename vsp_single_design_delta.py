import os
import openvsp as vsp
import pyvista as pv
import numpy as np
import pandas as pd
import csv
import glob
import plotly.graph_objects as go
import plotly.express as px
from dash import Dash, dcc, html, Input, Output
# import dash_bootstrap_components as dbc

vsp_exe = r"C:\Program Files\OpenVSP-3.47.0\vsp.exe"

wing_span_res = 20
wing_chord_res = 50
velocities = list(range(10, 25, 5)) # m/s
alphas = list(range(-5, 15, 3)) # degrees AoA

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
    Sref = 0.5 * (root_chord + root_chord * taper_ratio) * span
    bref = span
    cref = root_chord

    # Aero sweep
    csv_exists = False
    for v in velocities:
        print(f"\n=== Running VSP Aero Sweep at {v} m/s ===")
        stl_path, vsp3_path = generate_wing("wing")
        CL, CD, CDi, Cm, load_res = vsp_sweep(vsp3_path, v, Sref, bref, cref)
        cl_data = read_lift_distribution(load_res)
        first_aoa_key = list(cl_data.keys())[0]
        span_locations = cl_data[first_aoa_key]['SoverB']
        aero_headers = ["Velocity", "Alpha_deg", "CL", "CD", "Cm"] + [f"Cl_span_{loc:.4f}" for loc in span_locations[0]] + ["Oswald_efficiency"]

        # Combine sweep results with local spanwise results
        aero_results = []
        for i, alpha in enumerate(alphas):
            e = [compute_oswald(CL[i], CDi[i], Sref, bref)]
            base_row = [v, alpha, CL[i], CD[i], Cm[i]]
            closest_aoa_key = min(cl_data.keys(), key=lambda k: abs(k - alpha))
            spanwise_cls = cl_data[closest_aoa_key]['Cl'][0]
            full_row = base_row + list(spanwise_cls) + e
            aero_results.append(full_row)

        aero_filename = "aero_full.csv"
        with open(aero_filename, 'a', newline='') as f:
            writer = csv.writer(f)
            if not csv_exists:
                writer.writerow(aero_headers)
                csv_exists = True
            writer.writerows(aero_results)

        os.remove("w")
        for filename in glob.glob("wing*"):
            try:
                os.remove(filename)
            except OSError:
                pass

    # Stability sweep
    stab_csv_exists = False
    
    for v in velocities:
        stl_path, vsp3_path = generate_wing("wing")
        stab_results = vsp_stability(vsp3_path, v, Sref, bref, cref)
        stab_dict = read_stability(stab_results)
        sm = get_sm(stab_results)
        
        stability_filename = "stability.csv"
        stab_headers = ["Velocity"] + list(stab_dict.keys()) + ["StaticMargin"]
        stab_columns = [v] + list(stab_dict.values()) + [sm]
        
        with open(stability_filename, 'a', newline='') as f:
            writer = csv.writer(f)
            if not stab_csv_exists:
                writer.writerow(stab_headers)
                stab_csv_exists = True
            writer.writerow(stab_columns)

        os.remove("w")
        for filename in glob.glob("wing*"):
            try:
                os.remove(filename)
            except OSError:
                pass

def generate_wing(wing_name):
    tip_chord = root_chord * taper_ratio
    airfoil_fwd = airfoil_file.replace("\\", "/")

    vsp.VSPCheckSetup()
    vsp.ClearVSPModel()
    wing_id = vsp.AddGeom( "WING" )
    vsp.SetParmVal( wing_id, "TotalSpan",      "WingGeom", span)
    vsp.SetParmVal( wing_id, "Root_Chord",     "XSec_1",   root_chord)
    vsp.SetParmVal( wing_id, "Tip_Chord",      "XSec_1",   tip_chord)
    vsp.SetParmVal( wing_id, "Sweep",          "XSec_1",   sweep)
    vsp.SetParmVal( wing_id, "Dihedral",       "XSec_1",   dihedral)
    vsp.SetParmVal( wing_id, "Twist",          "XSec_1",   twist)
    vsp.SetParmVal( wing_id, "Twist_Location", "XSec_1",   1.0 )
    vsp.SetParmVal( wing_id, "SectTess_U",     "XSec_1",   wing_span_res)
    vsp.SetParmVal( wing_id, "Tess_W",         "Shape",    wing_chord_res)
    
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
    vsp.SetDoubleAnalysisInput(aero_analysis, "AlphaStart", [alphas[0]])
    vsp.SetIntAnalysisInput(aero_analysis, "AlphaNpts", [len(alphas)])
    vsp.SetDoubleAnalysisInput(aero_analysis, "AlphaEnd", [alphas[-1]])
    vsp.SetDoubleAnalysisInput(aero_analysis, "MachStart", [mach])
    vsp.SetIntAnalysisInput(aero_analysis, "MachNpts", [1])
    vsp.SetIntAnalysisInput(aero_analysis, "WakeNumIter", [6]) 
    vsp.SetDoubleAnalysisInput(aero_analysis, "Vinf", [v])
    vsp.SetDoubleAnalysisInput(aero_analysis, "Xcg", [x_cg])
    vsp.SetIntAnalysisInput(aero_analysis, "NCPU", [8])
    vsp.SetStringAnalysisInput(aero_analysis, "RedirectFile", f"{vsp3_path}_log.txt")
    rid = vsp.ExecAnalysis(aero_analysis)

    # Results
    polar_res = vsp.FindLatestResultsID("VSPAERO_Polar")
    load_res = vsp.FindLatestResultsID("VSPAERO_Load")
    cl = vsp.GetDoubleResults(polar_res, "CLtot")
    cd = vsp.GetDoubleResults(polar_res, "CDtot")
    cdi = vsp.GetDoubleResults(polar_res, "CDi")
    cm = vsp.GetDoubleResults(polar_res, "CMytot")
    return cl, cd, cdi, cm, load_res

def vsp_stability(vsp3_path, v, Sref, bref, cref):    
    # Load model
    vsp.ReadVSPFile(vsp3_path)
    geom_analysis = "VSPAEROComputeGeometry"
    vsp.SetAnalysisInputDefaults(geom_analysis)
    vsp.SetIntAnalysisInput(geom_analysis, "GeomSet", [vsp.SET_NONE])      
    vsp.SetIntAnalysisInput(geom_analysis, "ThinGeomSet", [vsp.SET_SHOWN])
    
    # Control surfaces
    cs_pitch_id = vsp.CreateVSPAEROControlSurfaceGroup() # Pitch
    vsp.SetVSPAEROControlGroupName("Pitch", cs_pitch_id)
    vsp.AddAllToVSPAEROControlSurfaceGroup(cs_pitch_id)
    cs_roll_id = vsp.CreateVSPAEROControlSurfaceGroup() # Roll
    vsp.SetVSPAEROControlGroupName("Roll", cs_roll_id)
    vsp.AddAllToVSPAEROControlSurfaceGroup(cs_roll_id)
    vsp.Update()
    group_pitch_str = f"ControlSurfaceGroup_{cs_pitch_id + 1}"
    container_id = vsp.FindContainer("VSPAEROSettings", 0)
    wing_id = vsp.FindGeoms()[0]
    cs_id = vsp.GetSubSurf(wing_id, 0)
    vsp.SetParmVal(vsp.FindParm(container_id, f"Surf_{cs_id}_1_Gain", group_pitch_str), -1)
    vsp.ExecAnalysis(geom_analysis)

    # Stability Sweep
    aero_analysis = "VSPAEROSweep"
    vsp.SetAnalysisInputDefaults(aero_analysis)
    vsp.SetDoubleAnalysisInput(aero_analysis, "Sref", [Sref])
    vsp.SetDoubleAnalysisInput(aero_analysis, "cref", [cref])
    vsp.SetDoubleAnalysisInput(aero_analysis, "bref", [bref])
    vsp.SetDoubleAnalysisInput(aero_analysis, "AlphaStart", [0.0])
    vsp.SetIntAnalysisInput(aero_analysis, "AlphaNpts", [1])
    vsp.SetDoubleAnalysisInput(aero_analysis, "AlphaEnd", [0.0])
    mach = v / 343.0
    vsp.SetDoubleAnalysisInput(aero_analysis, "MachStart", [mach])
    vsp.SetIntAnalysisInput(aero_analysis, "MachNpts", [1])
    vsp.SetIntAnalysisInput(aero_analysis, "WakeNumIter", [6]) 
    vsp.SetDoubleAnalysisInput(aero_analysis, "Vinf", [100.0])
    vsp.SetDoubleAnalysisInput(aero_analysis, "Xcg", [x_cg])
    vsp.SetIntAnalysisInput(aero_analysis, "UnsteadyType", [1])
    vsp.SetIntAnalysisInput(aero_analysis, "NCPU", [8])
    vsp.SetStringAnalysisInput(aero_analysis, "RedirectFile", f"{vsp3_path}_log.txt")
    vsp.ExecAnalysis(aero_analysis)
    
    stab_results = vsp.FindLatestResultsID("VSPAERO_Stab")
    return stab_results

def read_stability(stab_res):
    vsp_dict = {}
    
    # Stability Derivatives
    vsp_dict["Cm_alpha"] = vsp.GetDoubleResults(stab_res, "Alpha_CMm")[0] # Longitudinal Static Stability (C_m_alpha)
    vsp_dict["Cn_beta"]  = vsp.GetDoubleResults(stab_res, "Beta_CMn")[0] # Directional Static Stability (C_n_beta)
    vsp_dict["Cl_beta"]  = vsp.GetDoubleResults(stab_res, "Beta_CMl")[0] # Lateral Static Stability / Dihedral Effect (C_l_beta)
    vsp_dict["Cm_q"]     = vsp.GetDoubleResults(stab_res, "Pitch_Rate_CMm")[0] # Pitch Damping (C_m_q)
    vsp_dict["Cn_r"]     = vsp.GetDoubleResults(stab_res, "Yaw___Rate_CMn")[0] # Yaw Damping (C_n_r)
    vsp_dict["Cl_p"]     = vsp.GetDoubleResults(stab_res, "Roll__Rate_CMl")[0] # Roll Damping (C_l_p)
    vsp_dict["Cl_r"]     = vsp.GetDoubleResults(stab_res, "Yaw___Rate_CMl")[0] # Roll due to Yaw Rate (C_l_r) - Major Dutch Roll contributor
    vsp_dict["Cn_p"]     = vsp.GetDoubleResults(stab_res, "Roll__Rate_CMn")[0] # Yaw due to Roll Rate (C_n_p)

    # Control Derivatives (Safely fall back to 0.0 if geometry is missing)
    try:
        vsp_dict["Cm_de"] = vsp.GetDoubleResults(stab_res, "Pitch_CMm")[0] # Elevator Authority (C_m_de)
    except IndexError:
        vsp_dict["Cm_de"] = 0.0

    try:
        vsp_dict["Cl_da"] = vsp.GetDoubleResults(stab_res, "Roll_CMl")[0] # Aileron Authority (C_l_da)
    except IndexError:
        vsp_dict["Cl_da"] = 0.0

    try:
        vsp_dict["Cn_da"] = vsp.GetDoubleResults(stab_res, "Roll_CMn")[0] # Adverse Yaw (C_n_da)
    except IndexError:
        vsp_dict["Cn_da"] = 0.0

    return vsp_dict
    
def read_lift_distribution(load_result, target_vortex_sheet=1):
    data_by_aoa = {}
    current_aoa = vsp.GetDoubleResults(load_result, "Angle ")[0]
    if current_aoa not in data_by_aoa:
        data_by_aoa[current_aoa] = {'SoverB': [], 'Cl': []}
        soverb = vsp.GetDoubleResults(load_result, "SoverB")
        cl = vsp.GetDoubleResults(load_result, "cl*c/cref")
        data_by_aoa[current_aoa]['SoverB'].append(soverb)
        data_by_aoa[current_aoa]['Cl'].append(cl)         
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
    
def get_sm(stab_res):
    return vsp.GetDoubleResults(stab_res, "SM")[0]

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
        html.H1("Flying Wing Performance", className="mt-4 mb-4 text-center"),
        
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
                                value=0,
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
                # Row 1
                dbc.Row([
                    dbc.Col(dcc.Graph(id='stab-cl-beta'), width=4),
                    dbc.Col(dcc.Graph(id='stab-cn-beta'), width=4),
                    dbc.Col(dcc.Graph(id='stab-static-margin'), width=4),
                ], className="mt-3"),
                # Row 2
                dbc.Row([
                    dbc.Col(dcc.Graph(id='stab-cl-p'), width=4),
                    dbc.Col(dcc.Graph(id='stab-cn-p'), width=4),
                    dbc.Col(dcc.Graph(id='stab-cm-q'), width=4),
                ], className="mt-3"),
                # Row 3
                dbc.Row([
                    dbc.Col(dcc.Graph(id='stab-cl-r'), width=4),
                    dbc.Col(dcc.Graph(id='stab-cn-r'), width=4),
                    dbc.Col(dcc.Graph(id='stab-cm-de'), width=4),
                ], className="mt-3"),
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
            fig = px.line(df_aero, x=x_col, y=y_col, color='Velocity', title=title, markers=True, template=PLOT_TEMPLATE)
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
        [Output('stab-cl-beta', 'figure'), Output('stab-cn-beta', 'figure'), Output('stab-static-margin', 'figure'),
         Output('stab-cl-p', 'figure'), Output('stab-cn-p', 'figure'), Output('stab-cm-q', 'figure'),
         Output('stab-cl-r', 'figure'), Output('stab-cn-r', 'figure'), Output('stab-cm-de', 'figure')],
        [Input('alpha-dropdown', 'value')] # Using any input to trigger initial load
    )
    def update_stability_grid(_):
        def create_stab_fig(y_col, title, color="#3498DB"):
            fig = px.line(df_stab, x='Velocity', y=y_col, markers=True, template=PLOT_TEMPLATE, title=title)
            fig.update_traces(line_color=color)
            fig.update_layout(margin=dict(l=20, r=20, t=40, b=20))
            return fig

        return (
            create_stab_fig('Cl_beta', 'Cl_beta (Dihedral Effect)', '#E74C3C'),
            create_stab_fig('Cn_beta', 'Cn_beta (Directional Stability)', '#F1C40F'),
            create_stab_fig('StaticMargin', 'Static Margin', "#2ECC70"),
            create_stab_fig('Cl_p', 'Cl_p (Roll Damping)', "#2E65CC"),
            create_stab_fig('Cn_p', 'Cn_p (Cross Derivative)', "#F0F0F0"),
            create_stab_fig('Cm_q', 'Cm_q (Pitch Damping)', "#921CA1"),
            create_stab_fig('Cl_r', 'Cl_r (Cross Derivative)', "#7E510F"),
            create_stab_fig('Cn_r', 'Cn_r (Yaw Damping)', "#B95674"),
            create_stab_fig('Cm_de', 'Cm_de (Elevator Authority)', "#0DB9B9")
        )

    return app

if __name__ == '__main__':
    main()    
    # app = plot_dashboards()
    # app.run(debug=True)