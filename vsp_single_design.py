import os
import openvsp as vsp
import pyvista as pv
import os
import csv
import glob

wing_span_res = 20
wing_chord_res = 50
velocity = 10 # m/s
alphas = list(range(-15, 21)) # degrees AoA

airfoil_file = r"C:\Users\kprah\Desktop\Prahaas\WatArrow\CFD Automation\Airfoils\dae21.dat"

root_chord = 91.0
taper_ratio = 0.1
sweep = 55.0
dihedral = 5.0
twist = -15.0
span = 121.0
x_cg = 45.0
elevon_length = 0.1         # % of chord
elevon_start = 0.2          # % of wingspan
elevon_end = 0.8            # % of wingspan

def main():
    # Aero sweep
    stl_path, analysis_path = generate_wing("wing", span, root_chord, taper_ratio, sweep, dihedral, twist, airfoil_file, elevon_length, elevon_start, elevon_end)
    visualize_stl(stl_path)

    CL, CD, Cm = vsp_sweep(analysis_path, velocity, alphas, 0.5 * (root_chord + root_chord * taper_ratio) * span, span, root_chord)
    aero_results = zip(alphas, CL, CD, Cm)

    for filename in glob.glob(f"wing*"):
        try:
            os.remove(filename)
        except OSError:
            pass
    
    aero_filename = "cfd_sweep.csv"
    aero_headers = ["Alpha_deg", "CL", "CD", "Cm"]
    with open(aero_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(aero_headers)     # Write the title row
        writer.writerows(aero_results)    # Write all data rows

    # Stability Sweep
    stl_path, analysis_path = generate_wing("wing", span, root_chord, taper_ratio, sweep, dihedral, twist, airfoil_file, elevon_length, elevon_start, elevon_end)
    vsp_stability(analysis_path, velocity, [0.0], 0.5 * (root_chord + root_chord * taper_ratio) * span, span, root_chord)

    os.rename('wing.stab', 'STABILITY.txt')
    read_stability()
    os.remove("STABILITY.txt")

    for filename in glob.glob(f"wing*"):
        try:
            os.remove(filename)
        except OSError:
            pass

def generate_wing(wing_name, wingspan, root_chord, taper_ratio, sweep_angle, dihedral_angle, twist_angle, airfoil_file, elev_l, elev_s, elev_e):
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

    # Control surfaces
    cs_id = vsp.AddSubSurf(wing_id, vsp.SS_CONTROL)
    vsp.SetParmVal(wing_id, "SE_Const_Flag", "SS_Control_1", 1.0)
    vsp.SetParmVal(wing_id, "Length_C_Start", "SS_Control_1", elev_l)
    vsp.SetParmVal(wing_id, "EtaFlag", "SS_Control_1", 1.0)
    vsp.SetParmVal(wing_id, "EtaStart", "SS_Control_1", elev_s)
    vsp.SetParmVal(wing_id, "EtaEnd", "SS_Control_1", elev_e)

    # Finalize and export
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
    vsp.SetDoubleAnalysisInput(aero_analysis, "Vinf", [100.0])
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
    vsp.ReadVSPFile(fname_vspaerotests)
    
    # Geometry
    geom_analysis = "VSPAEROComputeGeometry"
    vsp.SetAnalysisInputDefaults(geom_analysis)
    
    vsp.SetIntAnalysisInput(geom_analysis, "GeomSet", [vsp.SET_NONE])      
    vsp.SetIntAnalysisInput(geom_analysis, "ThinGeomSet", [vsp.SET_SHOWN])
    
    # Control surfaces
    cs_pitch_id = vsp.CreateVSPAEROControlSurfaceGroup() # Pitch
    vsp.SetVSPAEROControlGroupName("Pitch", cs_pitch_id)
    cs_roll_id = vsp.CreateVSPAEROControlSurfaceGroup() # Roll
    vsp.SetVSPAEROControlGroupName("Roll", cs_roll_id)
    vsp.AddAllToVSPAEROControlSurfaceGroup(cs_pitch_id)
    vsp.AddAllToVSPAEROControlSurfaceGroup(cs_roll_id)

    container_id = vsp.FindContainer("VSPAEROSettings", 0)
    wing_id = vsp.FindGeoms()[0]
    cs_id = vsp.GetSubSurf(wing_id, 0)        
    vsp.SetParmVal(container_id, f"Surf_{cs_id}_0_Gain", "ControlSurfaceGroup_0", 1.0)
    vsp.SetParmVal(container_id, f"Surf_{cs_id}_1_Gain", "ControlSurfaceGroup_0", -1.0)
    vsp.SetParmVal(container_id, f"Surf_{cs_id}_0_Gain", "ControlSurfaceGroup_1", 1.0)
    vsp.SetParmVal(container_id, f"Surf_{cs_id}_1_Gain", "ControlSurfaceGroup_1", 1.0)
    
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
    vsp.SetDoubleAnalysisInput(aero_analysis, "Vinf", [100.0])
    vsp.SetDoubleAnalysisInput(aero_analysis, "Xcg", [x_cg])

    vsp.SetIntAnalysisInput(aero_analysis, "UnsteadyType", [1])

    print(f"--- Running Stability Sweep ({aero_analysis}) ---")
    vsp.ExecAnalysis(aero_analysis)

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
                    vsp_dict["CL_de"]   = float(parts[9])  # Extra lift from elevator
                elif coef == "CS":  # Side Force (Y)
                    vsp_dict["CY_beta"] = float(parts[3])
                    vsp_dict["CY_p"]    = float(parts[4])
                    vsp_dict["CY_r"]    = float(parts[6])
                elif coef == "CMl": # Roll Moment (l)
                    vsp_dict["Cl_beta"] = float(parts[3])
                    vsp_dict["Cl_p"]    = float(parts[4])
                    vsp_dict["Cl_r"]    = float(parts[6])
                    vsp_dict["Cl_da"]   = float(parts[10]) # Roll from Aileron
                elif coef == "CMm": # Pitch Moment (m)
                    vsp_dict["Cm_q"]    = float(parts[5])
                    vsp_dict["Cm_de"]   = float(parts[9])  # Pitch from Elevator
                elif coef == "CMn": # Yaw Moment (n)
                    vsp_dict["Cn_beta"] = float(parts[3])
                    vsp_dict["Cn_p"]    = float(parts[4])
                    vsp_dict["Cn_r"]    = float(parts[6])
                    vsp_dict["Cn_da"]   = float(parts[10]) # Adverse Yaw from Aileron
                
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