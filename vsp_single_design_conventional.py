import os
import openvsp as vsp
import pyvista as pv
import csv
import glob

wing_span_res = 20
wing_chord_res = 50
velocity = 10 # m/s
alphas = list(range(-15, 21)) # degrees AoA

airfoil_file = r"C:\Users\kprah\Desktop\Prahaas\WatArrow\CFD Automation\Airfoils\s1223.dat"

x_cg = 0.06

wing_params = {
    "span": 1.0,         # [m]
    "root_chord": 0.231, # [m]
    "taper": 1.0,        # [Ratio]
    "sweep": 0.0,        # [deg] Leading Edge Sweep
    "dihedral": 0.0,     # [deg]
    "twist": 0.0,        # [deg] Washout at tip
    "alpha": 0.0         # [deg]
}

htail_params = {
    "chord": 0.142,
    "l_H": 0.379,            # [mm] Tail Moment Arm (Distance from CG to Tail AC)
    "airfoil": "0012",   
    "span": 0.25,
    "alpha": 0.0 
}

vtail_params = {
    "chord": 0.142,
    "airfoil": "0012",   
    "span": 0.18,
    "taper": 0.75,
    "sweep": 20.0
}

def main():
    # Aero sweep
    # stl_path, analysis_path = generate_wing_and_tail("plane", wing_params, airfoil_file, htail_params)
    # visualize_stl(stl_path)

    # CL, CD, Cm = vsp_sweep(analysis_path, velocity, [0.0], wing_params["root_chord"] * wing_params["span"], wing_params["span"], wing_params["root_chord"])
    # aero_results = zip(alphas, CL, CD, Cm)

    # for filename in glob.glob(f"wing*"):
    #     try:
    #         os.remove(filename)
    #     except OSError:
    #         pass
    
    # aero_filename = "cfd_sweep.csv"
    # aero_headers = ["Alpha_deg", "CL", "CD", "Cm"]
    # with open(aero_filename, 'w', newline='') as f:
    #     writer = csv.writer(f)
    #     writer.writerow(aero_headers)     # Write the title row
    #     writer.writerows(aero_results)    # Write all data rows

    # Stability Sweep
    stl_path, analysis_path = generate_wing_and_tail("plane", wing_params, airfoil_file, htail_params)
    # visualize_stl(stl_path)
    vsp_stability(analysis_path, velocity, [0.0], wing_params["root_chord"] * wing_params["span"], wing_params["span"], wing_params["root_chord"])
    
    # os.rename('plane.stab', 'STABILITY.txt')
    # read_stability()
    # os.remove("STABILITY.text")

    # for filename in glob.glob(f"plane*"):
    #     try:
    #         os.remove(filename)
    #     except OSError:
    #         pass

def generate_wing_and_tail(wing_name, wing_params, airfoil_file, htail_params):
    vsp.VSPCheckSetup()
    vsp.ClearVSPModel()
    
    # Main wing
    wing_id = vsp.AddGeom("WING")
    vsp.SetGeomName(wing_id, "Main_Wing")
    tip_chord = wing_params["root_chord"] * wing_params["taper"]
    vsp.SetParmVal(wing_id, "TotalSpan", "WingGeom", wing_params["span"])
    vsp.SetParmVal(wing_id, "Root_Chord", "XSec_1", wing_params["root_chord"])
    vsp.SetParmVal(wing_id, "Tip_Chord", "XSec_1", tip_chord)
    vsp.SetParmVal(wing_id, "Sweep", "XSec_1", wing_params["sweep"])
    vsp.SetParmVal(wing_id, "Dihedral", "XSec_1", wing_params["dihedral"])
    vsp.SetParmVal(wing_id, "Twist", "XSec_1", wing_params["twist"])
    vsp.SetParmVal(wing_id, "Twist_Location", "XSec_1", 1.0) # Twist at tip
    vsp.SetParmVal(wing_id, "Y_Rel_Rotation", "XForm", wing_params["alpha"])
    vsp.SetParmVal(wing_id, "SectTess_U", "XSec_1", wing_span_res)
    vsp.SetParmVal(wing_id, "Tess_W", "Shape", wing_chord_res)
    
    for i in [0, 1]: # 0=Root, 1=Tip
        surf = vsp.GetXSecSurf(wing_id, i)
        vsp.ChangeXSecShape(surf, i, vsp.XS_FILE_AIRFOIL)
        xsec = vsp.GetXSec(surf, i)
        vsp.ReadFileAirfoil(xsec, airfoil_file)
        
    # Ailerons
    cs_a_id = vsp.AddSubSurf(wing_id, vsp.SS_CONTROL)
    vsp.SetParmVal(wing_id, "SE_Const_Flag", "SS_Control_1", 1.0)
    vsp.SetParmVal(wing_id, "Length_C_Start", "SS_Control_1", 0.2)
    vsp.SetParmVal(wing_id, "EtaFlag", "SS_Control_1", 1.0)
    vsp.SetParmVal(wing_id, "EtaStart", "SS_Control_1", 0.2)
    vsp.SetParmVal(wing_id, "EtaEnd", "SS_Control_1", 0.8)
        
    vsp.SetSetFlag(wing_id, 1, True)
    
    # Horizontal tail
    htail_id = vsp.AddGeom("WING")
    vsp.SetGeomName(htail_id, "Horizontal_Tail")
    vsp.SetParmVal(htail_id, "TotalSpan", "WingGeom", htail_params["span"])
    vsp.SetParmVal(htail_id, "Root_Chord", "XSec_1", htail_params["chord"])
    vsp.SetParmVal(htail_id, "Tip_Chord", "XSec_1", htail_params["chord"])
    vsp.SetParmVal(htail_id, "Sweep", "XSec_1", 0.0)
    vsp.SetParmVal(htail_id, "Dihedral", "XSec_1", 0.0)
    vsp.SetParmVal(htail_id, "Twist", "XSec_1", 0.0)
    vsp.SetParmVal(htail_id, "SectTess_U", "XSec_1", wing_span_res)
    vsp.SetParmVal(htail_id, "Tess_W", "Shape", wing_chord_res)
    vsp.SetParmVal(htail_id, "X_Rel_Location", "XForm", htail_params["l_H"])
    vsp.SetParmVal(htail_id, "Y_Rel_Rotation", "XForm", htail_params["alpha"])
    
    naca_code = htail_params["airfoil"]
    camber = int(naca_code[0]) / 100.0
    cam_loc = int(naca_code[1]) / 10.0
    thick = int(naca_code[2:]) / 100.0
    for i in [0, 1]:
        surf = vsp.GetXSecSurf(htail_id, i)
        # vsp.ChangeXSecShape(surf, i, vsp.XS_NACA_4_SERIES)
        xsec = vsp.GetXSec(surf, i)
        vsp.SetParmVal(htail_id, "Camber", f"XSecCurve_{i}", camber)
        vsp.SetParmVal(htail_id, "CamberLoc", f"XSecCurve_{i}", cam_loc)
        vsp.SetParmVal(htail_id, "ThickChord", f"XSecCurve_{i}", thick)
        
    # Elevator
    cs_e_id = vsp.AddSubSurf(htail_id, vsp.SS_CONTROL)
    vsp.SetParmVal(htail_id, "SE_Const_Flag", "SS_Control_1", 1.0)
    vsp.SetParmVal(htail_id, "Length_C_Start", "SS_Control_1", 0.25)
    vsp.SetParmVal(htail_id, "EtaFlag", "SS_Control_1", 1.0)
    vsp.SetParmVal(htail_id, "EtaStart", "SS_Control_1", 0.0)
    vsp.SetParmVal(htail_id, "EtaEnd", "SS_Control_1", 0.8)

    vsp.SetSetFlag(htail_id, 2, True)
    
    # Vertical tail
    vtail_id = vsp.AddGeom("WING")
    vsp.SetGeomName(vtail_id, "Vertical_Tail")
    vsp.SetParmVal(vtail_id, "Sym_Planar_Flag", "Sym", 0.0)
    vsp.SetParmVal(vtail_id, "TotalSpan", "WingGeom", vtail_params["span"])
    vsp.SetParmVal(vtail_id, "Root_Chord", "XSec_1", vtail_params["chord"])
    vsp.SetParmVal(vtail_id, "Tip_Chord", "XSec_1", vtail_params["chord"] * vtail_params["taper"])
    vsp.SetParmVal(vtail_id, "SectTess_U", "XSec_1", wing_span_res)
    vsp.SetParmVal(vtail_id, "Tess_W", "Shape", wing_chord_res)
    vsp.SetParmVal(vtail_id, "Sweep", "XSec_1", vtail_params["sweep"])
    vsp.SetParmVal(vtail_id, "Dihedral", "XSec_1", 0.0)
    vsp.SetParmVal(vtail_id, "Twist", "XSec_1", 0.0)
    vsp.SetParmVal(vtail_id, "X_Rel_Location", "XForm", htail_params["l_H"])
    vsp.SetParmVal(vtail_id, "X_Rel_Rotation", "XForm", 90.0)
    
    naca_code = vtail_params["airfoil"]
    camber = int(naca_code[0]) / 100.0
    cam_loc = int(naca_code[1]) / 10.0
    thick = int(naca_code[2:]) / 100.0
    for i in [0, 1]:
        surf = vsp.GetXSecSurf(vtail_id, i)
        # vsp.ChangeXSecShape(surf, i, vsp.XS_NACA_4_SERIES)
        xsec = vsp.GetXSec(surf, i)
        vsp.SetParmVal(htail_id, "Camber", f"XSecCurve_{i}", camber)
        vsp.SetParmVal(htail_id, "CamberLoc", f"XSecCurve_{i}", cam_loc)
        vsp.SetParmVal(htail_id, "ThickChord", f"XSecCurve_{i}", thick)
    
    # Rudder
    cs_r_id = vsp.AddSubSurf(vtail_id, vsp.SS_CONTROL)
    vsp.SetParmVal(vtail_id, "SE_Const_Flag", "SS_Control_1", 1.0)
    vsp.SetParmVal(vtail_id, "Length_C_Start", "SS_Control_1", 0.35)
    vsp.SetParmVal(vtail_id, "EtaFlag", "SS_Control_1", 1.0)
    vsp.SetParmVal(vtail_id, "EtaStart", "SS_Control_1", 0.2)
    vsp.SetParmVal(vtail_id, "EtaEnd", "SS_Control_1", 0.8)

    vsp.SetSetFlag(vtail_id, 3, True)

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

def vsp_sweep(fname_vspaerotests, vin, alphas, Sref, bref, cref):
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
    vsp.SetDoubleAnalysisInput(aero_analysis, "AlphaStart", [alphas[0]])
    vsp.SetIntAnalysisInput(aero_analysis, "AlphaNpts", [len(alphas)])
    vsp.SetDoubleAnalysisInput(aero_analysis, "AlphaEnd", [alphas[-1]])    
    mach = vin / 343.0
    vsp.SetDoubleAnalysisInput(aero_analysis, "MachStart", [mach])
    vsp.SetIntAnalysisInput(aero_analysis, "MachNpts", [1])
    vsp.SetIntAnalysisInput(aero_analysis, "WakeNumIter", [15]) 
    vsp.SetDoubleAnalysisInput(aero_analysis, "Vinf", [vin])
    vsp.SetDoubleAnalysisInput(aero_analysis, "Xcg", [x_cg])

    print(f"--- Running Aero Sweep ({aero_analysis}) ---")
    rid = vsp.ExecAnalysis(aero_analysis)

    # Results
    polar_res = vsp.FindLatestResultsID("VSPAERO_Polar")
    cl = vsp.GetDoubleResults(polar_res, "CLtot")
    cd = vsp.GetDoubleResults(polar_res, "CDtot")
    cm = vsp.GetDoubleResults(polar_res, "CMytot")
    return cl, cd, cm

def vsp_stability(fname_vspaerotests, vin, alphas, Sref, bref, cref):
    # Load model
    vsp.ClearVSPModel()
    vsp.ReadVSPFile(fname_vspaerotests)
    
    # Geometry
    geom_analysis = "VSPAEROComputeGeometry"
    vsp.SetAnalysisInputDefaults(geom_analysis)
    
    vsp.SetIntAnalysisInput(geom_analysis, "GeomSet", [vsp.SET_NONE])      
    vsp.SetIntAnalysisInput(geom_analysis, "ThinGeomSet", [vsp.SET_ALL])
    
    # Control surfaces
    vsp.AutoGroupVSPAEROControlSurfaces()    
    vsp.Update()

    container_id = vsp.FindContainer("VSPAEROSettings", 0)
    elev_id = vsp.FindGeoms()[1]
    cs_id = vsp.GetSubSurf(elev_id, 0)
    vsp.SetParmVal(vsp.FindParm(container_id, f"Surf_{cs_id}_1_Gain", "Horizontal_Tail_SS_CONT_0"), -1) # Since elevators are symmetrical
    
    vsp.Update()
    print(f"--- Running Meshing ({geom_analysis}) ---")
    vsp.ExecAnalysis(geom_analysis)

    aero_analysis = "VSPAEROSweep"
    vsp.SetAnalysisInputDefaults(aero_analysis)
    
    # Reference dimensions
    vsp.SetDoubleAnalysisInput(aero_analysis, "Sref", [Sref])
    vsp.SetDoubleAnalysisInput(aero_analysis, "cref", [cref])
    vsp.SetDoubleAnalysisInput(aero_analysis, "bref", [bref])
    
    # Flight conditions    
    vsp.SetDoubleAnalysisInput(aero_analysis, "AlphaStart", [alphas[0]])
    vsp.SetIntAnalysisInput(aero_analysis, "AlphaNpts", [len(alphas)])
    vsp.SetDoubleAnalysisInput(aero_analysis, "AlphaEnd", [alphas[-1]])    
    mach = vin / 343.0
    vsp.SetDoubleAnalysisInput(aero_analysis, "MachStart", [mach])
    vsp.SetIntAnalysisInput(aero_analysis, "MachNpts", [1])
    vsp.SetIntAnalysisInput(aero_analysis, "WakeNumIter", [15]) 
    vsp.SetDoubleAnalysisInput(aero_analysis, "Vinf", [vin])
    vsp.SetDoubleAnalysisInput(aero_analysis, "Xcg", [x_cg])

    vsp.SetIntAnalysisInput(aero_analysis, "UnsteadyType", [1])
    
    vsp.SetDoubleAnalysisInput(aero_analysis, "Rho", [1.225])
    
    vsp.SetStringAnalysisInput(aero_analysis, "RedirectFile", ["log.txt"])
    vsp.SetIntAnalysisInput(aero_analysis, "NCPU", [8])

    print(f"--- Running Stability Sweep ({aero_analysis}) ---")
    vsp.ExecAnalysis(aero_analysis)
    
    stab_results = vsp.FindLatestResultsID("VSPAERO_Stab")
    vsp.PrintResults(stab_results)

def read_stability(input_file="STABILITY.txt", output_file="vsp_derivatives.csv"):
    target_aoa = "0.0000000" 
    current_aoa = None
    capture_derivatives = False
    vsp_dict = {}

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue

            if line.startswith("AoA_"):
                current_aoa = line.split()[1]
                continue

            if "Coef" in line and "Total" in line:
                capture_derivatives = True
                continue

            if capture_derivatives and current_aoa == target_aoa:
                parts = line.split()
                if not parts: continue
                
                coef = parts[0]
                
                if coef == "CL": # Lift 
                    vsp_dict["CL_de"]   = float(parts[10])  # Extra lift from elevator
                elif coef == "CS":  # Side Force (Y)
                    vsp_dict["CY_beta"] = float(parts[3])
                    vsp_dict["CY_p"]    = float(parts[4])
                    vsp_dict["CY_r"]    = float(parts[6])
                elif coef == "CMl": # Roll Moment (l)
                    vsp_dict["Cl_beta"] = float(parts[3])
                    vsp_dict["Cl_p"]    = float(parts[4])
                    vsp_dict["Cl_r"]    = float(parts[6])
                    vsp_dict["Cl_da"]   = float(parts[9]) # Roll from Aileron
                elif coef == "CMm": # Pitch Moment (m)
                    vsp_dict["Cm_q"]    = float(parts[5])
                    vsp_dict["Cm_de"]   = float(parts[10])  # Pitch from Elevator
                elif coef == "CMn": # Yaw Moment (n)
                    vsp_dict["Cn_beta"] = float(parts[3])
                    vsp_dict["Cn_p"]    = float(parts[4])
                    vsp_dict["Cn_r"]    = float(parts[6])
                    vsp_dict["Cn_da"]   = float(parts[9]) # Adverse Yaw from Aileron
                
                elif "Result" in line or line.startswith("****"):
                    capture_derivatives = False

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Derivative", "Value"])
        for key, val in vsp_dict.items():
            writer.writerow([key, val])

    print(f"Successfully generated {output_file} formatted for the 6-DOF simulator.")

if __name__ == "__main__":
    main()