import os
import openvsp as vsp # type: ignore
import pyvista as pv
import numpy as np
import pandas as pd
import csv
import glob
import plotly.express as px
from dash import Dash, dcc, html

vsp_exe = r"C:\Program Files\OpenVSP-3.47.0\vsp.exe"

wing_span_res = 10
wing_chord_res = 25
velocities = list(range(10, 50, 5)) # m/s
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
    Sref = 0.5 * (root_chord + root_chord * taper_ratio) * span
    bref = span
    cref = root_chord

    # Aero sweep
    csv_exists = False
    for v in velocities:
        print(f"\n=== Running VSP Aero Sweep at {v} m/s ===")
        stl_path, vsp3_path = generate_wing("wing")
        CL, CD, CDi, Cm = vsp_sweep(vsp3_path, v, Sref, bref, cref)
        lift = [0.5 * 1.225 * (v ** 2) * Sref * cl_val for cl_val in CL]
        drag = [0.5 * 1.225 * (v ** 2) * Sref * cd_val for cd_val in CD]
        aero_headers = ["Velocity", "Alpha_deg", "CL", "CD", "Cm", "Lift", "Drag", "Oswald_efficiency"]

        # Combine sweep results with local spanwise results
        aero_results = []
        for i, alpha in enumerate(alphas):
            row = [v, alpha, CL[i], CD[i], Cm[i], lift[i], drag[i], compute_oswald(CL[i], CDi[i], Sref, bref)]
            aero_results.append(row)

        aero_filename = "aero_full.csv"
        with open(aero_filename, 'a', newline='') as f:
            writer = csv.writer(f)
            if not csv_exists:
                writer.writerow(aero_headers)
                csv_exists = True
            writer.writerows(aero_results)

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
        stab_dict = read_stability()
        
        stability_filename = "stability.csv"
        stab_headers = ["Velocity"] + list(stab_dict.keys())
        stab_columns = [v] + list(stab_dict.values())
        
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
    
    # Control surfaces
    elevon_id = vsp.AddSubSurf(wing_id, vsp.SS_CONTROL)
    vsp.SetSubSurfName(elevon_id, "Elevons")
    vsp.SetParmVal(wing_id, "EtaFlag", "SS_Control_1", 1.0)
    vsp.SetParmVal(wing_id, "EtaStart", "SS_Control_1", elevon_start)
    vsp.SetParmVal(wing_id, "EtaEnd", "SS_Control_1", elevon_end)
    vsp.SetParmVal(wing_id, "SE_Const_Flag", "SS_Control_1", 1.0)
    vsp.SetParmVal(wing_id, "Length_C_Start", "SS_Control_1", elevon_length)
    cs_pitch_id = vsp.CreateVSPAEROControlSurfaceGroup() # Pitch
    vsp.SetVSPAEROControlGroupName("Pitch", cs_pitch_id)
    vsp.AddAllToVSPAEROControlSurfaceGroup(cs_pitch_id)
    cs_roll_id = vsp.CreateVSPAEROControlSurfaceGroup() # Roll
    vsp.SetVSPAEROControlGroupName("Roll", cs_roll_id)
    vsp.AddAllToVSPAEROControlSurfaceGroup(cs_roll_id)
    vsp.Update()
    container_id = vsp.FindContainer("VSPAEROSettings", 0)
    vsp.SetParmVal(vsp.FindParm(container_id, f"Surf_{elevon_id}_0_Gain", "ControlSurfaceGroup_0"), 1.0)
    vsp.SetParmVal(vsp.FindParm(container_id, f"Surf_{elevon_id}_1_Gain", "ControlSurfaceGroup_0"), -1.0)
    vsp.SetParmVal(vsp.FindParm(container_id, f"Surf_{elevon_id}_0_Gain", "ControlSurfaceGroup_1"), 1.0)
    vsp.SetParmVal(vsp.FindParm(container_id, f"Surf_{elevon_id}_1_Gain", "ControlSurfaceGroup_1"), 1.0)

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
    vsp.SetIntAnalysisInput(geom_analysis, "ThinGeomSet", [vsp.SET_SHOWN])
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
    vsp.SetStringAnalysisInput(aero_analysis, "RedirectFile", [f"{vsp3_path}_log.txt"])
    rid = vsp.ExecAnalysis(aero_analysis)

    # Results
    polar_res = vsp.FindLatestResultsID("VSPAERO_Polar")
    cl = vsp.GetDoubleResults(polar_res, "CLtot")
    cd = vsp.GetDoubleResults(polar_res, "CDtot")
    cdi = vsp.GetDoubleResults(polar_res, "CDi")
    cm = vsp.GetDoubleResults(polar_res, "CMytot")
    return cl, cd, cdi, cm

def vsp_stability(vsp3_path, v, Sref, bref, cref):
    # Load model
    vsp.ReadVSPFile(vsp3_path)
    geom_analysis = "VSPAEROComputeGeometry"
    vsp.SetAnalysisInputDefaults(geom_analysis)   
    vsp.SetIntAnalysisInput(geom_analysis, "GeomSet", [vsp.SET_NONE])      
    vsp.SetIntAnalysisInput(geom_analysis, "ThinGeomSet", [vsp.SET_SHOWN])
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
    vsp.SetStringAnalysisInput(aero_analysis, "RedirectFile", [f"{vsp3_path}_log.txt"])
    vsp.ExecAnalysis(aero_analysis)

def read_stability(stab_file_path="wing.stab"):
    vsp_dict = {
        "Cm_alpha": 0.0, "Cn_beta": 0.0, "Cl_beta": 0.0,
        "Cm_q": 0.0, "Cn_r": 0.0, "Cl_p": 0.0,
        "Cl_r": 0.0, "Cn_p": 0.0, "Static Margin": 0.0,
        "Cm_de": 0.0, "Cl_da": 0.0, "Cn_da": 0.0
    }

    with open(stab_file_path, 'r') as f:
        lines = f.readlines()

    header_idx = -1
    for i, line in enumerate(lines):
        if line.startswith("Coef") and "Total" in line and "Alpha" in line:
            header_idx = i
            break

    if header_idx != -1:
        headers = lines[header_idx].split()
        
        col = {
            "Alpha": headers.index("Alpha"),
            "Beta": headers.index("Beta"),
            "p": headers.index("p"),
            "q": headers.index("q"),
            "r": headers.index("r"),
            "ConGrp_1": headers.index("ConGrp_1"),
            "ConGrp_2": headers.index("ConGrp_2")
        }

        for line in lines[header_idx+1:]:
            parts = line.split()
            if not parts or parts[0].startswith("#"):
                continue
            
            coef = parts[0]
            try:
                if coef == "CMl":
                    vsp_dict["Cl_beta"] = float(parts[col["Beta"]])
                    vsp_dict["Cl_p"]    = float(parts[col["p"]])
                    vsp_dict["Cl_r"]    = float(parts[col["r"]])
                    vsp_dict["Cl_da"]   = float(parts[col["ConGrp_2"]]) # Roll authority
                elif coef == "CMm":
                    vsp_dict["Cm_alpha"] = float(parts[col["Alpha"]])
                    vsp_dict["Cm_q"]     = float(parts[col["q"]])
                    vsp_dict["Cm_de"]    = float(parts[col["ConGrp_1"]]) # Pitch authority
                elif coef == "CMn":
                    vsp_dict["Cn_beta"]  = float(parts[col["Beta"]])
                    vsp_dict["Cn_p"]     = float(parts[col["p"]])
                    vsp_dict["Cn_r"]     = float(parts[col["r"]])
                    vsp_dict["Cn_da"]    = float(parts[col["ConGrp_2"]]) # Adverse Yaw
            except IndexError:
                pass

    for line in reversed(lines):
        if line.startswith("SM"):
            parts = line.split()
            if len(parts) >= 2:
                vsp_dict["Static Margin"] = float(parts[1])
            break

    return vsp_dict
    
def read_lift_distribution(load_result):
    data_by_aoa = {}
    for current_aoa in vsp.GetDoubleResults(load_result, "FC_AoA_"):
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
    
def plot_dashboards(sweep_csv="aero_full.csv", stab_csv="stability.csv"):
    df_aero = pd.read_csv(sweep_csv)
    df_stab = pd.read_csv(stab_csv)
    df_aero['L_D'] = df_aero['CL'] / df_aero['CD']

    # Initialize Dash App
    app = Dash(__name__)
    PLOT_TEMPLATE = "plotly_dark"

    # --- Figure Generation Helpers ---
    def create_aero_fig(y_col, title, x_col='Alpha_deg'):
        fig = px.line(df_aero, x=x_col, y=y_col, color='Velocity', title=title, markers=True, template=PLOT_TEMPLATE)
        fig.update_layout(margin=dict(l=20, r=20, t=40, b=20), paper_bgcolor='#222222', plot_bgcolor='#222222')
        return fig

    def create_stab_fig(y_col, title, color="#3498DB"):
        # Fallback if a column hasn't been generated yet
        if y_col not in df_stab.columns:
            fig = px.line(title=f"{title} (No Data)", template=PLOT_TEMPLATE)
        else:
            fig = px.line(df_stab, x='Velocity', y=y_col, markers=True, template=PLOT_TEMPLATE, title=title)
            fig.update_traces(line_color=color)
        
        fig.update_layout(margin=dict(l=20, r=20, t=40, b=20), paper_bgcolor='#222222', plot_bgcolor='#222222')
        return fig

    # --- CSS Styles ---
    grid_style = {
        'display': 'grid',
        'gridTemplateColumns': 'repeat(3, 1fr)', # 3 columns of equal width
        'gap': '20px',
        'padding': '20px'
    }
    
    tab_style = {'backgroundColor': '#333333', 'color': 'white', 'border': 'none'}
    tab_selected_style = {'backgroundColor': '#555555', 'color': 'white', 'fontWeight': 'bold', 'border': 'none'}

    # --- App Layout ---
    app.layout = html.Div(style={'backgroundColor': '#111111', 'minHeight': '100vh', 'color': 'white', 'fontFamily': 'Arial, sans-serif'}, children=[
        html.H1("Flying Wing Performance", style={'textAlign': 'center', 'paddingTop': '20px', 'margin': '0'}),
        
        dcc.Tabs(style={'padding': '20px'}, children=[
            
            # --- TAB 1: AERODYNAMICS ---
            dcc.Tab(label='Aerodynamics', style=tab_style, selected_style=tab_selected_style, children=[
                html.Div(style=grid_style, children=[
                    dcc.Graph(figure=create_aero_fig('CL', 'CL vs Alpha')),
                    dcc.Graph(figure=create_aero_fig('CD', 'CD vs Alpha')),
                    dcc.Graph(figure=create_aero_fig('Cm', 'Cm vs Alpha')),
                    dcc.Graph(figure=create_aero_fig('Lift', 'Lift vs Alpha')),
                    dcc.Graph(figure=create_aero_fig('Drag', 'Drag vs Alpha')),
                    dcc.Graph(figure=create_aero_fig('CL', 'Drag Polar (CL vs CD)', x_col='CD')),
                    dcc.Graph(figure=create_aero_fig('L_D', 'L/D vs Alpha')),
                    dcc.Graph(figure=create_aero_fig('Oswald_efficiency', 'Oswald Efficiency vs Alpha'))
                ])
            ]),

            # --- TAB 2: STABILITY ---
            dcc.Tab(label='Stability', style=tab_style, selected_style=tab_selected_style, children=[
                html.Div(style=grid_style, children=[
                    dcc.Graph(figure=create_stab_fig('Cl_beta', 'Cl_beta (Dihedral Effect)', '#E74C3C')),
                    dcc.Graph(figure=create_stab_fig('Cn_beta', 'Cn_beta (Directional Stability)', '#F1C40F')),
                    dcc.Graph(figure=create_stab_fig('Static Margin', 'Static Margin', "#2ECC70")),
                    dcc.Graph(figure=create_stab_fig('Cl_p', 'Cl_p (Roll Damping)', "#2E65CC")),
                    dcc.Graph(figure=create_stab_fig('Cn_p', 'Cn_p (Cross Derivative)', "#F0F0F0")),
                    dcc.Graph(figure=create_stab_fig('Cm_q', 'Cm_q (Pitch Damping)', "#921CA1")),
                    dcc.Graph(figure=create_stab_fig('Cl_r', 'Cl_r (Cross Derivative)', "#7E510F")),
                    dcc.Graph(figure=create_stab_fig('Cn_r', 'Cn_r (Yaw Damping)', "#B95674")),
                    dcc.Graph(figure=create_stab_fig('Cm_de', 'Cm_de (Elevator Authority)', "#0DB9B9"))
                ])
            ])
        ])
    ])

    return app

if __name__ == '__main__':
    # main()
    app = plot_dashboards()
    app.run(debug=True)