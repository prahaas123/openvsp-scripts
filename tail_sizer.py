import os
import openvsp as vsp
import pyvista as pv
import math
import uuid
import glob
import os

wing_span_res = 20
wing_chord_res = 50
velocity = 10 # m/s
alpha = 0 # degrees AoA

airfoil_file = r"C:\Users\kprah\Desktop\Prahaas\WatArrow\CFD Automation\Airfoils\s1223.dat"

moment_tolerance = 0.05
tail_sizing_iterations = 10

wing_params = {
    "span": 1.0,          # [m]
    "root_chord": 0.231,  # [m]
    "taper": 1.0,         # [Ratio]
    "sweep": 0.0,         # [deg] Leading Edge Sweep
    "dihedral": 0.0,      # [deg]
    "twist": 0.0,         # [deg] Washout at tip
    "alpha": 0.0          # [deg]
}

htail_params = {
    "V_H": 0.1,           # Horizontal Tail Volume Coefficient (Typical: 0.4 - 0.6)
    "l_H": 0.379,         # [cm] Tail Moment Arm (Distance from CG to Tail AC)
    "airfoil": "0012",   
    "aspect_ratio": 3.5  
}

vtail_params = {
    "V_V": 0.04,         # Vertical Tail Volume Coefficient (Typical: 0.03 - 0.06)
    "airfoil": "0012",
    "aspect_ratio": 1.5
}

def main():
    # Generate just main wing
    name = f"wing_{uuid.uuid4().hex[:8]}"
    generate_wing(name, wing_params, airfoil_file)
    plain_moment = get_moment(name, velocity, alpha, 0.25 * wing_params["root_chord"], wing_params["span"] * wing_params["root_chord"], wing_params["root_chord"])
    print(f"Pitching Moment Coefficient (Cmy) at alpha={alpha} deg: {plain_moment:.6f} (without a tail)")
    
    # Generate initial tail geometry
    tail_name = f"plane_{uuid.uuid4().hex[:8]}"
    S_tail = (htail_params["V_H"] * wing_params["span"] * wing_params["span"] * wing_params["root_chord"]) / htail_params["l_H"]
    b_tail = math.sqrt(S_tail * htail_params["aspect_ratio"])    
    generate_wing_and_htail(tail_name, wing_params, airfoil_file, htail_params, b_tail, 0)
    moment = get_moment(tail_name, velocity, alpha, 0.25 * wing_params["root_chord"], wing_params["span"] * wing_params["root_chord"], wing_params["root_chord"])
    
    # Tail sizing loop
    i_old = 0.0
    m_old = moment
    success = False
    print(f"Iter 0: Tail Alpha = {i_old:.2f} deg -> CMy = {m_old:.6f}")
    
    if abs(m_old) < moment_tolerance:
        print("Aircraft is naturally trimmed!")
        success = True
    else:
        i_new = -2.0 # Initial guess
        tail_name = f"plane_{uuid.uuid4().hex[:8]}"
        generate_wing_and_htail(tail_name, wing_params, airfoil_file, htail_params, b_tail, i_new)
        m_new = get_moment(tail_name, velocity, alpha, 0.25 * wing_params["root_chord"], wing_params["span"] * wing_params["root_chord"], wing_params["root_chord"])
        print(f"Iter 1: Tail Alpha = {i_new:.2f} deg -> CMy = {m_new:.6f}")
        
        for iteration in range(2, tail_sizing_iterations + 1):
            if abs(m_new) < moment_tolerance:
                print(f"SUCCESS: Trimmed at Tail Alpha = {i_new:.4f} degrees")
                success = True
                break
            if abs(m_new - m_old) < 1e-9:
                print("Slope is zero! Cannot converge further.")
                break
            i_next = i_new - m_new * (i_new - i_old) / (m_new - m_old)
            i_next = max(min(i_next, 15.0), -15.0)
            i_old = i_new
            m_old = m_new
            i_new = i_next
            
            tail_name = f"plane_{uuid.uuid4().hex[:8]}"
            generate_wing_and_htail(tail_name, wing_params, airfoil_file, htail_params, b_tail, i_new)
            m_new = get_moment(tail_name, velocity, alpha, 0.25 * wing_params["root_chord"], wing_params["span"] * wing_params["root_chord"], wing_params["root_chord"])
            print(f"Iter {iteration}: Tail Alpha = {i_new:.2f} deg -> CMy = {m_new:.6f}")
        else:
            print("WARNING: Max iterations reached without full convergence.")
    
    if success:
        print(f"Final Tail Alpha for Trim: {i_new:.4f} degrees with CMy = {m_new:.6f}")
        print(f"Horizontal tail dimensions: Chord: {b_tail / htail_params["aspect_ratio"]:.3f}m, Span: {b_tail:.3f}m")
        generate_wing_and_htail("plane_final", wing_params, airfoil_file, htail_params, b_tail, i_new)
        visualize_stl("plane_final.stl")

def generate_wing(wing_name, wing_params, airfoil_file):
    vsp.VSPCheckSetup()
    vsp.ClearVSPModel()

    # Create the wing
    wing_id = vsp.AddGeom("WING")
    
    # Setting parameters
    tip_chord = wing_params["root_chord"] * wing_params["taper"]
    
    vsp.SetParmVal(wing_id, "TotalSpan", "WingGeom", wing_params["span"])
    vsp.SetParmVal(wing_id, "Root_Chord", "XSec_1", wing_params["root_chord"])
    vsp.SetParmVal(wing_id, "Tip_Chord", "XSec_1", tip_chord)
    vsp.SetParmVal(wing_id, "Sweep", "XSec_1", wing_params["sweep"])
    vsp.SetParmVal(wing_id, "Dihedral", "XSec_1", wing_params["dihedral"])
    vsp.SetParmVal(wing_id, "Twist", "XSec_1", wing_params["twist"])
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

def generate_wing_and_htail(wing_name, wing_params, airfoil_file, htail_params, htail_b, htail_alpha):
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
        
    vsp.SetSetFlag(wing_id, 1, True)

    # Horizontal tail
    htail_id = vsp.AddGeom("WING")
    vsp.SetGeomName(htail_id, "Horizontal_Tail")
    tail_chord = htail_b / htail_params["aspect_ratio"]
    vsp.SetParmVal(htail_id, "TotalSpan", "WingGeom", htail_b)
    vsp.SetParmVal(htail_id, "Root_Chord", "XSec_1", tail_chord)
    vsp.SetParmVal(htail_id, "Tip_Chord", "XSec_1", tail_chord)
    vsp.SetParmVal(htail_id, "Sweep", "XSec_1", 0.0)
    vsp.SetParmVal(htail_id, "Dihedral", "XSec_1", 0.0)
    vsp.SetParmVal(htail_id, "Twist", "XSec_1", 0.0)
    vsp.SetParmVal(htail_id, "SectTess_U", "XSec_1", wing_span_res)
    vsp.SetParmVal(htail_id, "Tess_W", "Shape", wing_chord_res)
    vsp.SetParmVal(htail_id, "X_Rel_Location", "XForm", htail_params["l_H"])
    vsp.SetParmVal(htail_id, "Y_Rel_Rotation", "XForm", htail_alpha)

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

    vsp.SetSetFlag(htail_id, 1, True) # Add to Set 1

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

def get_moment(fname_vspaerotests, velocity, alpha, x_cg, Sref, cref):
    vsp.ReadVSPFile(f"{fname_vspaerotests}.vsp3")
    
    # Geometry
    geom_analysis = "VSPAEROComputeGeometry"
    vsp.SetAnalysisInputDefaults(geom_analysis)
    
    vsp.SetIntAnalysisInput(geom_analysis, "GeomSet", [vsp.SET_NONE])      
    vsp.SetIntAnalysisInput(geom_analysis, "ThinGeomSet", [vsp.SET_SHOWN])

    vsp.ExecAnalysis(geom_analysis)

    aero_analysis = "VSPAEROSweep"
    vsp.SetAnalysisInputDefaults(aero_analysis)
    
    # Reference dimensions
    vsp.SetDoubleAnalysisInput(aero_analysis, "Sref", [Sref])
    vsp.SetDoubleAnalysisInput(aero_analysis, "cref", [cref])
    vsp.SetDoubleAnalysisInput(aero_analysis, "Xcg", [x_cg])
    
    # Flight conditions    
    vsp.SetDoubleAnalysisInput(aero_analysis, "AlphaStart", [alpha])
    vsp.SetIntAnalysisInput(aero_analysis, "AlphaNpts", [1])    
    mach = velocity / 343.0
    vsp.SetDoubleAnalysisInput(aero_analysis, "MachStart", [mach])
    vsp.SetIntAnalysisInput(aero_analysis, "MachNpts", [1])
    vsp.SetIntAnalysisInput(aero_analysis, "WakeNumIter", [15]) 
    vsp.SetDoubleAnalysisInput(aero_analysis, "Vinf", [100.0])
    
    vsp.SetIntAnalysisInput(aero_analysis, "NCPU", [8])
    rid = vsp.ExecAnalysis(aero_analysis)

    # Results
    polar_res = vsp.FindLatestResultsID("VSPAERO_Polar")        
    moment = vsp.GetDoubleResults(polar_res, "CMytot")[0]
    
    for filename in glob.glob(f"{fname_vspaerotests}*"):
            try:
                os.remove(filename)
            except OSError:
                pass
            
    return moment
    
if __name__ == "__main__":
    main()
    # S_tail = (htail_params["V_H"] * wing_params["span"] * wing_params["span"] * wing_params["root_chord"]) / htail_params["l_H"]
    # b_tail = math.sqrt(S_tail * htail_params["aspect_ratio"])
    # print(f"Initial tail area: {S_tail:.2f} cm^2, Tail span: {b_tail:.2f} cm, Initial tail chord: {S_tail / b_tail:.2f} cm")