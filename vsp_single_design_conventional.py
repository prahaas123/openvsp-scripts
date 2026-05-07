import os
import pyvista as pv
import csv
import glob
import numpy as np
import openvsp as vsp # type: ignore
import shutil

wing_span_res = 10
wing_chord_res = 25
velocities = list(range(10, 50, 5)) # m/s
alphas = list(range(-5, 15)) # degrees AoA

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
        aero_headers = ["Velocity", "Alpha_deg", "CL", "CD", "Cm", "Oswald_efficiency"]

        aero_results = []
        for i, alpha in enumerate(alphas):
            e = compute_oswald(CL[i], CDi[i], Sref, bref)
            row = [v, alpha, CL[i], CD[i], Cm[i], e]
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

    # Stability sweep
    stab_csv_exists = False
    
    for v in velocities:
        vsp3_path = shutil.copy(PLANE, "plane.vsp3")
        vsp_stability(vsp3_path, v, Sref, bref, cref)
        stab_dict = read_stability("plane.stab")
        
        stability_filename = "stability.csv"
        stab_headers = ["Velocity"] + list(stab_dict.keys())
        stab_columns = [v] + list(stab_dict.values())
        
        with open(stability_filename, 'a', newline='') as f:
            writer = csv.writer(f)
            if not stab_csv_exists:
                writer.writerow(stab_headers)
                stab_csv_exists = True
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
    
    for i in [0, 1]:
        surf = vsp.GetXSecSurf(wid, i)
        vsp.ChangeXSecShape(surf, 0, vsp.XS_FILE_AIRFOIL)
        vsp.ReadFileAirfoil(vsp.GetXSec(surf, 0), airfoil_fwd)
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
            if tokens[0] == 'Coef':
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
 
    vsp_dict = {
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

if __name__ == "__main__":
    main()