import os
import csv
import glob
import numpy as np
import openvsp as vsp # type: ignore
import shutil
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

wing_span_res = 10
wing_chord_res = 25
velocities = list(range(10, 50, 5)) # m/s
alphas = list(range(-5, 15)) # degrees AoA

stability_velocity = 20.0 # m/s
airfoil_file = r"Airfoils\goe322.dat"
x_cg = 0.065

wing_params = {
    "span": 1.19,          # [m]
    "root_chord": 0.225,   # [m]
    "taper": 1.0,         # [Ratio]
    "sweep": 0.0,         # [deg] Leading Edge Sweep
    "dihedral": 0.0,      # [deg]
    "twist": 0.0,         # [deg] Washout at tip
    "alpha": 3.0          # [deg]
}

htail_params = {
    "chord": 0.128,
    "l_H": 0.65,        # [m] Tail Moment Arm (Distance from CG to Tail AC)
    "airfoil": "0012",
    "span": 0.383,
    "alpha": 3.35
}

vtail_params = {
    "chord": 0.128,
    "airfoil": "0012",
    "span": 0.18,
    "taper": 0.7,
    "sweep": 15.0
}

def main():
    bref = wing_params["span"]
    cref = wing_params["root_chord"]
    Sref = 0.5 * (cref + (cref * wing_params["taper"])) * bref
    
    _, PLANE = generate_wing_and_tail("Conventional")
    
    try:
        os.remove("aero_full.csv")
        os.remove("stability.csv")
    except OSError:
        pass

    # Aero sweep
    csv_exists = False
    for v in velocities:
        print(f"\n=== Running VSP Aero Sweep at {v} m/s ===")
        vsp3_path = shutil.copy(PLANE, "plane.vsp3")
        CL, CD, CDi, Cm = vsp_sweep(vsp3_path, v, Sref, bref, cref)
        aero_headers = ["Velocity", "Alpha_deg", "CL", "CD", "Cm", "Lift", "Drag", "Oswald_efficiency"]

        aero_results = []
        for i, alpha in enumerate(alphas):
            e = compute_oswald(CL[i], CDi[i], Sref, bref)
            lift = 0.5 * 1.225 * (v ** 2) * Sref * CL[i]
            drag = 0.5 * 1.225 * (v ** 2) * Sref * CD[i]
            row = [v, alpha, CL[i], CD[i], Cm[i], lift, drag, e]
            aero_results.append(row)

        # Write to CSV
        aero_filename = "aero_full.csv"
        with open(aero_filename, 'a', newline='') as f:
            writer = csv.writer(f)
            if not csv_exists:
                writer.writerow(aero_headers)
                csv_exists = True
            writer.writerows(aero_results)

        # File Cleanup
        for filename in glob.glob("plane*"):
            try:
                os.remove(filename)
            except OSError:
                pass

    # Stability
    v = stability_velocity
    print(f"\n=== Running VSP Stability at {v} m/s (cruise) ===")
    vsp3_path = shutil.copy(PLANE, "plane.vsp3")
    vsp_stability(vsp3_path, v, Sref, bref, cref)
    stab_dict = read_stability("plane.stab")

    stability_filename = "stability.csv"
    stab_headers = ["Velocity"] + list(stab_dict.keys())
    stab_columns = [v] + list(stab_dict.values())

    with open(stability_filename, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(stab_headers)
        writer.writerow(stab_columns)

    for filename in glob.glob("plane*"):
        try:
            os.remove(filename)
        except OSError:
            pass
            
    for filename in glob.glob("Conventional*"):
        try:
            os.remove(filename)
        except OSError:
            pass

    # Launch dashboard
    try:
        plot_dashboard("aero_full.csv", "stability.csv")
    except Exception as exc:
        print(f"[dashboard] skipped: {exc}")

def generate_wing_and_tail(plane_name):
    airfoil_fwd  = airfoil_file.replace("\\", "/")
    tip_chord    = wing_params["root_chord"] * wing_params["taper"]
    vtail_tip    = vtail_params["chord"] * vtail_params["taper"]

    def naca4(code):
        return int(code[0]) / 100.0, int(code[1]) / 10.0, int(code[2:]) / 100.0
    h_camber, h_cam_loc, h_thick = naca4(htail_params["airfoil"])
    v_camber, v_cam_loc, v_thick = naca4(vtail_params["airfoil"])

    # Main Wing
    wid = vsp.AddGeom("WING", "")
    vsp.SetGeomName(wid, "MainWing")
    vsp.SetParmVal(wid, "TotalSpan", "WingGeom", wing_params["span"])
    vsp.SetParmVal(wid, "Root_Chord", "XSec_1", wing_params["root_chord"])
    vsp.SetParmVal(wid, "Tip_Chord", "XSec_1", tip_chord)
    vsp.SetParmVal(wid, "Sweep", "XSec_1", wing_params["sweep"])
    vsp.SetParmVal(wid, "Dihedral", "XSec_1", wing_params["dihedral"])
    vsp.SetParmVal(wid, "Twist", "XSec_1", wing_params["twist"])
    vsp.SetParmVal(wid, "Y_Rel_Rotation", "XForm", wing_params["alpha"])
    vsp.SetParmVal(wid, "SectTess_U", "XSec_1", float(wing_span_res))
    vsp.SetParmVal(wid, "Tess_W", "Shape", float(wing_chord_res))
    
    surf = vsp.GetXSecSurf(wid, 0)
    for i in [0, 1]:
        vsp.ChangeXSecShape(surf, i, vsp.XS_FILE_AIRFOIL)
        vsp.ReadFileAirfoil(vsp.GetXSec(surf, i), airfoil_file)
    vsp.SetSetFlag(wid, 1, True)

    # Horizontal Tail
    hid = vsp.AddGeom("WING", "")
    vsp.SetGeomName(hid, "HorizontalTail")
    vsp.SetParmVal(hid, "TotalSpan", "WingGeom", htail_params["span"])
    vsp.SetParmVal(hid, "Root_Chord", "XSec_1", htail_params["chord"])
    vsp.SetParmVal(hid, "Tip_Chord", "XSec_1", htail_params["chord"])
    vsp.SetParmVal(hid, "Sweep", "XSec_1", 0.0)
    vsp.SetParmVal(hid, "X_Rel_Location", "XForm", htail_params["l_H"])
    vsp.SetParmVal(hid, "Y_Rel_Rotation", "XForm", htail_params["alpha"])
    vsp.SetParmVal(hid, "Camber", "XSecCurve_0", h_camber)
    vsp.SetParmVal(hid, "CamberLoc", "XSecCurve_0", h_cam_loc)
    vsp.SetParmVal(hid, "ThickChord", "XSecCurve_0", h_thick)
    vsp.SetSetFlag(hid, 1, True)

    # Vertical Tail
    vid = vsp.AddGeom("WING", "")
    vsp.SetGeomName(vid, "VerticalTail")
    vsp.SetParmVal(vid, "Sym_Planar_Flag", "Sym", 0.0)
    vsp.SetParmVal(vid, "TotalSpan", "WingGeom", vtail_params["span"])
    vsp.SetParmVal(vid, "Root_Chord", "XSec_1", vtail_params["chord"])
    vsp.SetParmVal(vid, "Tip_Chord", "XSec_1", vtail_tip)
    vsp.SetParmVal(vid, "Sweep", "XSec_1", 0.0)
    vsp.SetParmVal(vid, "X_Rel_Location", "XForm", htail_params["l_H"])
    vsp.SetParmVal(vid, "X_Rel_Rotation", "XForm", 90.0)
    vsp.SetParmVal(vid, "Camber", "XSecCurve_0", v_camber)
    vsp.SetParmVal(vid, "CamberLoc", "XSecCurve_0", v_cam_loc)
    vsp.SetParmVal(vid, "ThickChord", "XSecCurve_0", v_thick)
    vsp.SetSetFlag(vid, 1, True)

    vsp.Update()
    vsp3_path = f"{plane_name}.vsp3"
    stl_path = f"{plane_name}.stl"
    vsp.WriteVSPFile(vsp3_path)
    vsp.ExportFile(stl_path, 0, vsp.EXPORT_STL)
    return stl_path, vsp3_path

def vsp_sweep(vsp3_path, velocity, Sref, bref, cref):
    mach = velocity / 343.0

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
    vsp.SetDoubleAnalysisInput(aero_analysis, "Vinf", [velocity])
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
    vsp.ClearVSPModel()
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
    vsp.SetDoubleAnalysisInput(aero_analysis, "Vinf", [v])
    vsp.SetDoubleAnalysisInput(aero_analysis, "Xcg", [x_cg])
    vsp.SetIntAnalysisInput(aero_analysis, "UnsteadyType", [1])
    vsp.SetIntAnalysisInput(aero_analysis, "NCPU", [8])
    vsp.SetStringAnalysisInput(aero_analysis, "RedirectFile", [f"{vsp3_path}_log.txt"])
    vsp.ExecAnalysis(aero_analysis)

def parse_polar(polar_path):
    CL, CD, CDi, Cm = [], [], [], []
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
                col_cdi = tokens.index('CDi')
                col_cm = tokens.index('CMytot')
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

def read_stability(stab_path): 
    col_alpha = col_beta = col_p = col_q = col_r = None
    col_aileron = col_elevator = col_rudder = None
    rows = {}
 
    with open(stab_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            tokens = stripped.split()
            if tokens[0] == 'Coef' and col_alpha is None:
                col_alpha    = tokens.index('Alpha')
                col_beta     = tokens.index('Beta')
                col_p        = tokens.index('p')
                col_q        = tokens.index('q')
                col_r        = tokens.index('r')
                col_aileron  = tokens.index('ConGrp_1') if 'ConGrp_1' in tokens else None  # MainWing_SS_CONT_0
                col_elevator = tokens.index('ConGrp_2') if 'ConGrp_2' in tokens else None  # HorizontalTail_SS_CONT_0
                col_rudder   = tokens.index('ConGrp_3') if 'ConGrp_3' in tokens else None  # VerticalTail_SS_CONT_0
                continue
            if col_alpha is None:
                continue
            coef = tokens[0]
            if coef in ('CL', 'CS', 'CMl', 'CMm', 'CMn'):
                try:
                    rows[coef] = tokens
                except (ValueError, IndexError):
                    pass
 
    def get(coef, col):
        if coef in rows and col is not None:
            try:
                return float(rows[coef][col])
            except (ValueError, IndexError):
                pass
        return float('nan')
 
    cm_alpha = get('CMm', col_alpha)           # CMm wrt Alpha
    cl_alpha = get('CL',  col_alpha)           # CL  wrt Alpha
    with np.errstate(divide='ignore', invalid='ignore'):
        sm = float(-cm_alpha / cl_alpha) if (np.isfinite(cl_alpha) and cl_alpha != 0.0) else float('nan')

    vsp_dict = {
        'Cm_alpha': cm_alpha,                  # CMm wrt Alpha (longitudinal static stability)
        'CL_alpha': cl_alpha,                  # CL  wrt Alpha
        'Static_Margin': sm,                   # -Cm_alpha / CL_alpha  (>0 = statically stable)
        'CL_de':   get('CL',  col_elevator),   # CL  wrt elevator
        'CY_beta': get('CS',  col_beta),       # CS  wrt Beta
        'CY_p':    get('CS',  col_p),          # CS  wrt p
        'CY_r':    get('CS',  col_r),          # CS  wrt r
        'Cl_beta': get('CMl', col_beta),       # CMl wrt Beta
        'Cl_p':    get('CMl', col_p),          # CMl wrt p
        'Cl_r':    get('CMl', col_r),          # CMl wrt r
        'Cl_da':   get('CMl', col_aileron),    # CMl wrt aileron
        'Cm_q':    get('CMm', col_q),          # CMm wrt q
        'Cm_de':   get('CMm', col_elevator),   # CMm wrt elevator
        'Cn_beta': get('CMn', col_beta),       # CMn wrt Beta
        'Cn_p':    get('CMn', col_p),          # CMn wrt p
        'Cn_r':    get('CMn', col_r),          # CMn wrt r
        'Cn_da':   get('CMn', col_aileron),    # CMn wrt aileron (adverse yaw)
    }
 
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

def plot_dashboard(sweep_csv="aero_full.csv", stab_csv="stability.csv"):
    pio.renderers.default = "browser"

    df = pd.read_csv(sweep_csv)
    df["L_D"] = df["CL"] / df["CD"]
    stab = pd.read_csv(stab_csv).iloc[0].to_dict()

    velocities_sorted = sorted(df["Velocity"].unique())
    scale = pio.templates
    n = max(len(velocities_sorted) - 1, 1)
    vel_colors = {
        v: pv_color for v, pv_color in zip(
            velocities_sorted,
            [f"hsl({int(210 + 130 * i / n)},70%,60%)" for i in range(len(velocities_sorted))]
        )
    }

    specs = [
        [{"type": "xy"}, {"type": "xy"}, {"type": "xy"}, {"type": "xy"}],   # aero 1-4
        [{"type": "xy"}, {"type": "xy"}, {"type": "xy"}, {"type": "xy"}],   # aero 5-8
        [{"type": "xy", "colspan": 4}, None, None, None],                   # stability cards
        [{"type": "table", "colspan": 4}, None, None, None],                # values table
    ]
    titles = [
        "CL vs \u03b1", "CD vs \u03b1", "Cm vs \u03b1", "Drag polar (CL vs CD)",
        "L/D vs \u03b1", "Lift vs \u03b1 (N)", "Drag vs \u03b1 (N)", "Oswald e vs \u03b1",
        "Cruise stability check (\u03b1 = 0\u00b0)",
        "Stability derivatives",
    ]
    fig = make_subplots(
        rows=4, cols=4, specs=specs, subplot_titles=titles,
        row_heights=[0.20, 0.20, 0.24, 0.36],
        vertical_spacing=0.06, horizontal_spacing=0.06,
    )

    # Aero graphs
    def add_curves(xcol, ycol, row, col, show_legend):
        for v in velocities_sorted:
            sub = df[df["Velocity"] == v].sort_values("Alpha_deg")
            fig.add_trace(go.Scatter(
                x=sub[xcol], y=sub[ycol], mode="lines+markers",
                name=f"{v:.0f} m/s", legendgroup=f"{v}", showlegend=show_legend,
                line=dict(color=vel_colors[v], width=2), marker=dict(size=4),
            ), row=row, col=col)

    add_curves("Alpha_deg", "CL", 1, 1, True)
    add_curves("Alpha_deg", "CD", 1, 2, False)
    add_curves("Alpha_deg", "Cm", 1, 3, False)
    add_curves("CD", "CL", 1, 4, False)          # drag polar
    add_curves("Alpha_deg", "L_D", 2, 1, False)
    add_curves("Alpha_deg", "Lift", 2, 2, False)
    add_curves("Alpha_deg", "Drag", 2, 3, False)
    add_curves("Alpha_deg", "Oswald_efficiency", 2, 4, False)

    _amin, _amax = df["Alpha_deg"].min(), df["Alpha_deg"].max()
    fig.add_trace(go.Scatter(
        x=[_amin, _amax], y=[0, 0], mode="lines", hoverinfo="skip", showlegend=False,
        line=dict(color="#888", dash="dot", width=1),
    ), row=1, col=3)

    # axis labels
    fig.update_xaxes(title_text="\u03b1 (deg)", row=1, col=1); fig.update_yaxes(title_text="CL", row=1, col=1)
    fig.update_xaxes(title_text="\u03b1 (deg)", row=1, col=2); fig.update_yaxes(title_text="CD", row=1, col=2)
    fig.update_xaxes(title_text="\u03b1 (deg)", row=1, col=3); fig.update_yaxes(title_text="Cm", row=1, col=3)
    fig.update_xaxes(title_text="CD", row=1, col=4);         fig.update_yaxes(title_text="CL", row=1, col=4)
    fig.update_xaxes(title_text="\u03b1 (deg)", row=2, col=1); fig.update_yaxes(title_text="L/D", row=2, col=1)
    fig.update_xaxes(title_text="\u03b1 (deg)", row=2, col=2); fig.update_yaxes(title_text="Lift (N)", row=2, col=2)
    fig.update_xaxes(title_text="\u03b1 (deg)", row=2, col=3); fig.update_yaxes(title_text="Drag (N)", row=2, col=3)
    fig.update_xaxes(title_text="\u03b1 (deg)", row=2, col=4); fig.update_yaxes(title_text="e", row=2, col=4)

    # Stability scores
    GREEN, RED, GREY = "#2ECC71", "#E74C3C", "#7F8C8D"
    GREEN_BG, RED_BG, GREY_BG = "rgba(46,204,113,0.10)", "rgba(231,76,60,0.10)", "rgba(127,140,141,0.08)"
    cards = [
        ("Cm_alpha",      "Cm_\u03b1",        "target < 0",      lambda x: x < 0),
        ("Static_Margin", "Static Margin",    "target 0.05\u20130.15", lambda x: x > 0),
        ("Cn_beta",       "Cn_\u03b2",        "target > 0",      lambda x: x > 0),
        ("Cl_beta",       "Cl_\u03b2",        "target < 0",      lambda x: x < 0),
        ("Cm_q",          "Cm_q",             "target < 0",      lambda x: x < 0),
        ("Cl_p",          "Cl_p",             "target < 0",      lambda x: x < 0),
        ("Cn_r",          "Cn_r",             "target < 0",      lambda x: x < 0),
        ("CY_beta",       "CY_\u03b2",        "target < 0",      lambda x: x < 0),
    ]

    fig.add_trace(go.Scatter(x=[0, 4], y=[0, 2], mode="markers", marker=dict(opacity=0), hoverinfo="skip", showlegend=False), row=3, col=1)
    fig.update_xaxes(visible=False, range=[0, 4], row=3, col=1)
    fig.update_yaxes(visible=False, range=[0, 2], row=3, col=1)
    for idx, (key, label, target, test) in enumerate(cards):
        cx = idx % 4
        top_band = 2 - (idx // 4)
        x0, x1 = cx + 0.04, cx + 0.96
        y0, y1 = (top_band - 1) + 0.10, top_band - 0.10
        xc = (x0 + x1) / 2.0
        val = float(stab.get(key, float("nan")))
        if not np.isfinite(val):
            edge, fillc, valtxt = GREY, GREY_BG, "N/A"
        else:
            ok = test(val)
            edge, fillc = (GREEN, GREEN_BG) if ok else (RED, RED_BG)
            valtxt = f"{val:.3f}"
        fig.add_shape(type="rect", x0=x0, y0=y0, x1=x1, y1=y1, line=dict(color=edge, width=2), fillcolor=fillc, layer="below", row=3, col=1)
        fig.add_annotation(x=xc, y=y1 - 0.16, text=label, showarrow=False, font=dict(size=15, color="#CCC"), row=3, col=1)
        fig.add_annotation(x=xc, y=(y0 + y1) / 2.0 - 0.02, text=valtxt, showarrow=False, font=dict(size=26, color=edge), row=3, col=1)
        fig.add_annotation(x=xc, y=y0 + 0.14, text=target, showarrow=False, font=dict(size=11, color="#888"), row=3, col=1)

    # Values table
    ref = {
        "Cm_alpha":      ("< 0",            lambda x: x < 0),
        "CL_alpha":      ("> 0 (~4\u20136 /rad)", lambda x: x > 0),
        "Static_Margin": ("0.05 \u2013 0.15", lambda x: x > 0),
        "CY_beta":       ("< 0",            lambda x: x < 0),
        "CY_p":          ("\u2248 0 (small)", None),
        "CY_r":          ("> 0",            lambda x: x > 0),
        "Cl_beta":       ("< 0 (mild)",     lambda x: x < 0),
        "Cl_p":          ("< 0",            lambda x: x < 0),
        "Cl_r":          ("> 0",            lambda x: x > 0),
        "Cm_q":          ("< 0 (strong)",   lambda x: x < 0),
        "Cn_beta":       ("> 0",            lambda x: x > 0),
        "Cn_p":          ("\u2248 0 (small)", None),
        "Cn_r":          ("< 0",            lambda x: x < 0),
        "CL_de":         ("\u2014",         None),
        "Cm_de":         ("\u2014",         None),
        "Cl_da":         ("\u2014",         None),
        "Cn_da":         ("\u2014",         None),
    }
    order = ["Cm_alpha", "CL_alpha", "Static_Margin", "CY_beta", "CY_p", "CY_r",
             "Cl_beta", "Cl_p", "Cl_r", "Cm_q", "Cn_beta", "Cn_p", "Cn_r",
             "CL_de", "Cm_de", "Cl_da", "Cn_da"]
    NEUTRAL, PASS_BG, FAIL_BG, NA_BG = "#242424", "#183a28", "#3a1e1e", "#333"
    names, vals, recs, valcol = [], [], [], []
    for key in order:
        v = float(stab.get(key, float("nan")))
        rec_text, test = ref.get(key, ("", None))
        names.append(key)
        recs.append(rec_text)
        if not np.isfinite(v):
            vals.append("N/A"); valcol.append(NA_BG)
        else:
            vals.append(f"{v:.4f}")
            valcol.append(NEUTRAL if test is None else (PASS_BG if test(v) else FAIL_BG))
    colcol = [NEUTRAL] * len(order)
    fig.add_trace(go.Table(
        columnwidth=[1.0, 1.0, 1.4],
        header=dict(values=["Derivative", "Value", "Recommended"], fill_color="#444", font=dict(color="white", size=13), align="left"),
        cells=dict(values=[names, vals, recs], fill_color=[colcol, valcol, colcol], font=dict(color="#EEE", size=12), align="left", height=24),
    ), row=4, col=1)

    v_cruise = stab.get("Velocity", "?")
    fig.update_layout(
        template="plotly_dark",
        title=dict(text=f"Aero Sweep & Stability "
                        f"(V = {v_cruise} m/s, \u03b1 = 0\u00b0)",
                   x=0.5, font=dict(size=22)),
        height=1900, width=1500,
        paper_bgcolor="#111", plot_bgcolor="#1b1b1b",
        legend=dict(title="Velocity", orientation="v", x=1.01, y=1.0),
        margin=dict(l=60, r=60, t=90, b=40),
    )
    fig.show()
    return fig

if __name__ == "__main__":
    main()