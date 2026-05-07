import os
import openvsp as vsp # type: ignore
import math
import uuid
import glob
import pyvista as pv

wing_span_res = 20
wing_chord_res = 50
velocity = 10 # m/s
alpha = 0 # degrees AoA
SM = 0.10  # Desired Static Margin

airfoil_file = r"Airfoils\goe322.dat"

moment_tolerance = 0.05
tail_sizing_iterations = 10

wing_params = {
    "span": 1.19,
    "root_chord": 0.225,
    "taper": 1.0,
    "sweep": 0.0,
    "dihedral": 0.0,
    "twist": 0.0,
    "alpha": 3.0
}

htail_params = {
    "V_H": 0.10,
    "l_H": 0.65,
    "airfoil": "0012",
    "aspect_ratio": 3.0
}

vtail_params = {
    "V_V": 0.05,
    "airfoil": "0012",
    "aspect_ratio": 1.5
} 

def main():
    vsp.VSPCheckSetup()
    Sref = wing_params["span"] * wing_params["root_chord"]
    cref = wing_params["root_chord"]
    x_cg = 0.25 * wing_params["root_chord"]

    # Generate just main wing to get initial moment
    name = f"wing_{uuid.uuid4().hex[:8]}"
    vsp3_path = generate_wing(name)
    plain_moment = get_moment(vsp3_path, x_cg, Sref, cref)
    print(f"Pitching Moment Coefficient (Cmy) at alpha={alpha} deg: {plain_moment:.6f} (without a tail)")

    # Generate initial tail geometry
    S_tail = (htail_params["V_H"] * wing_params["span"]**2 * wing_params["root_chord"]) / htail_params["l_H"]
    b_tail = math.sqrt(S_tail * htail_params["aspect_ratio"])
    vtail_height = math.sqrt((vtail_params["V_V"] * wing_params["span"]**2 * wing_params["root_chord"]) / htail_params["l_H"])

    tail_name = f"plane_{uuid.uuid4().hex[:8]}"
    new_x_cg = calc_cg(S_tail)
    _, vsp3_path = generate_wing_and_tail(tail_name, b_tail, 0.0, vtail_height)
    moment = get_moment(vsp3_path, new_x_cg, Sref, cref)

    # Tail incidence sizing loop (secant method)
    i_old = 0.0
    m_old = moment
    success = False
    print(f"Iter 0: Tail Alpha = {i_old:.2f} deg -> CMy = {m_old:.6f}")

    if abs(m_old) < moment_tolerance:
        print("Aircraft is naturally trimmed!")
        i_new = i_old
        m_new = m_old
        success = True
    else:
        i_new = -2.0  # Initial guess
        tail_name = f"plane_{uuid.uuid4().hex[:8]}"
        _, vsp3_path = generate_wing_and_tail(tail_name, b_tail, i_new, vtail_height)
        m_new = get_moment(vsp3_path, new_x_cg, Sref, cref)
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
            i_old, m_old = i_new, m_new
            i_new = i_next

            tail_name = f"plane_{uuid.uuid4().hex[:8]}"
            _, vsp3_path = generate_wing_and_tail(tail_name, b_tail, i_new, vtail_height)
            m_new = get_moment(vsp3_path, x_cg, Sref, cref)
            print(f"Iter {iteration}: Tail Alpha = {i_new:.2f} deg -> CMy = {m_new:.6f}")
        else:
            print("WARNING: Max iterations reached without full convergence.")

    if success:
        tail_chord = b_tail / htail_params["aspect_ratio"]
        print(f"Final Tail Alpha for Trim: {i_new:.4f} degrees with CMy = {m_new:.6f}")
        print(f"Horizontal tail dimensions: Chord: {tail_chord:.3f}m, Span: {b_tail:.3f}m")
        print(f"Final CG location: {new_x_cg:.3f}m from LE, for a Static Margin of {SM*100:.1f}%")
        stl_final, _ = generate_wing_and_tail("plane_final", b_tail, i_new, vtail_height)
        visualize_stl(stl_final)

def generate_wing(wing_name):
    vsp.ClearVSPModel()
    wing_id = vsp.AddGeom("WING")
    tip_chord = wing_params["root_chord"] * wing_params["taper"]

    vsp.SetParmVal(wing_id, "TotalSpan", "WingGeom", wing_params["span"])
    vsp.SetParmVal(wing_id, "Root_Chord", "XSec_1", wing_params["root_chord"])
    vsp.SetParmVal(wing_id, "Tip_Chord", "XSec_1", tip_chord)
    vsp.SetParmVal(wing_id, "Sweep", "XSec_1", wing_params["sweep"])
    vsp.SetParmVal(wing_id, "Dihedral", "XSec_1", wing_params["dihedral"])
    vsp.SetParmVal(wing_id, "Twist", "XSec_1", wing_params["twist"])
    vsp.SetParmVal(wing_id, "Twist_Location", "XSec_1", 1.0)
    vsp.SetParmVal(wing_id, "SectTess_U", "XSec_1", float(wing_span_res))
    vsp.SetParmVal(wing_id, "Tess_W", "Shape", float(wing_chord_res))

    for i in [0, 1]:
        surf = vsp.GetXSecSurf(wing_id, i)
        vsp.ChangeXSecShape(surf, 0, vsp.XS_FILE_AIRFOIL)
        xsec = vsp.GetXSec(surf, 0)
        vsp.ReadFileAirfoil(xsec, airfoil_file)

    vsp.SetSetFlag(wing_id, 1, True)
    vsp.Update()
    vsp3_path = f"{wing_name}.vsp3"
    vsp.WriteVSPFile(vsp3_path)
    return vsp3_path

def generate_wing_and_tail(plane_name, htail_b, htail_alpha, vtail_height):
    vsp.ClearVSPModel()
    tip_chord = wing_params["root_chord"] * wing_params["taper"]
    tail_chord = htail_b / htail_params["aspect_ratio"]

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
        vsp.ReadFileAirfoil(vsp.GetXSec(surf, 0), airfoil_file)
    vsp.SetSetFlag(wid, 1, True)

    # Horizontal Tail
    hid = vsp.AddGeom("WING", "")
    vsp.SetGeomName(hid, "HorizontalTail")
    vsp.SetParmVal(hid, "TotalSpan", "WingGeom", htail_b)
    vsp.SetParmVal(hid, "Root_Chord", "XSec_1", tail_chord)
    vsp.SetParmVal(hid, "Tip_Chord", "XSec_1", tail_chord)
    vsp.SetParmVal(hid, "Sweep", "XSec_1", 0.0)
    vsp.SetParmVal(hid, "X_Rel_Location", "XForm", htail_params["l_H"])
    vsp.SetParmVal(hid, "Y_Rel_Rotation", "XForm", htail_alpha)
    vsp.SetParmVal(hid, "Camber", "XSecCurve_0", h_camber)
    vsp.SetParmVal(hid, "ThickChord", "XSecCurve_0", h_thick)
    vsp.SetSetFlag(hid, 1, True)

    # Vertical Tail
    vid = vsp.AddGeom("WING", "")
    vsp.SetGeomName(vid, "VerticalTail")
    vsp.SetParmVal(vid, "Sym_Planar_Flag", "Sym", 0.0)
    vsp.SetParmVal(vid, "TotalSpan", "WingGeom", vtail_height)
    vsp.SetParmVal(vid, "Root_Chord", "XSec_1", tail_chord)
    vsp.SetParmVal(vid, "Tip_Chord", "XSec_1", tail_chord)
    vsp.SetParmVal(vid, "Sweep", "XSec_1", 0.0)
    vsp.SetParmVal(vid, "X_Rel_Location", "XForm", htail_params["l_H"])
    vsp.SetParmVal(vid, "X_Rel_Rotation", "XForm", 90.0)
    vsp.SetParmVal(vid, "Camber", "XSecCurve_0", v_camber)
    vsp.SetParmVal(vid, "ThickChord", "XSecCurve_0", v_thick)
    vsp.SetSetFlag(vid, 1, True)

    vsp.Update()
    vsp3_path = f"{plane_name}.vsp3"
    stl_path = f"{plane_name}.stl"
    vsp.WriteVSPFile(vsp3_path)
    vsp.ExportFile(stl_path, 0, vsp.EXPORT_STL)
    return stl_path, vsp3_path

def get_moment(vsp3_path, x_cg, Sref, cref):
    vsp.ClearVSPModel()
    vsp.ReadVSPFile(vsp3_path)
    mach = velocity / 343.0

    # Geometry Compute
    vsp.SetAnalysisInputDefaults("VSPAEROComputeGeometry")
    vsp.SetIntAnalysisInput("VSPAEROComputeGeometry", "GeomSet", [vsp.SET_NONE])
    vsp.SetIntAnalysisInput("VSPAEROComputeGeometry", "ThinGeomSet", [vsp.SET_ALL])
    vsp.ExecAnalysis("VSPAEROComputeGeometry")

    # Aero Sweep
    vsp.SetAnalysisInputDefaults("VSPAEROSweep")
    vsp.SetDoubleAnalysisInput("VSPAEROSweep", "Sref", [Sref])
    vsp.SetDoubleAnalysisInput("VSPAEROSweep", "cref", [cref])
    vsp.SetDoubleAnalysisInput("VSPAEROSweep", "bref", [wing_params["span"]])
    vsp.SetDoubleAnalysisInput("VSPAEROSweep", "Xcg", [x_cg])
    vsp.SetDoubleAnalysisInput("VSPAEROSweep", "AlphaStart", [float(alpha)])
    vsp.SetDoubleAnalysisInput("VSPAEROSweep", "AlphaEnd", [float(alpha)])
    vsp.SetIntAnalysisInput("VSPAEROSweep", "AlphaNpts", [1])
    vsp.SetDoubleAnalysisInput("VSPAEROSweep", "MachStart", [mach])
    vsp.SetDoubleAnalysisInput("VSPAEROSweep", "Vinf", [velocity])
    vsp.SetIntAnalysisInput("VSPAEROSweep", "WakeNumIter", [8])
    vsp.SetStringAnalysisInput("VSPAEROSweep", "RedirectFile", [f"{vsp3_path}_log.txt"])
    
    vsp.ExecAnalysis("VSPAEROSweep")

    # Extract Results
    res_id = vsp.FindLatestResultsID("VSPAERO_Polar")
    cm = vsp.GetDoubleResults(res_id, "CMytot")[0]
    
    # Cleanup files
    base = os.path.splitext(vsp3_path)[0]
    for f in glob.glob(f"{base}*"):
        try: os.remove(f)
        except: pass
        
    return cm

def calc_cg(S_tail, SM=SM):
    b, c_r, taper = wing_params["span"], wing_params["root_chord"], wing_params["taper"]
    l_H, AR_t = htail_params["l_H"], htail_params["aspect_ratio"]
    S_w = 0.5 * (c_r + taper * c_r) * b
    MAC = (2/3) * c_r * (1 + taper + taper**2) / (1 + taper)
    AR_w = b**2 / S_w
    
    # Simple estimate for lift curves
    a_w = (2 * math.pi * AR_w) / (2 + math.sqrt(4 + AR_w**2))
    a_t = (2 * math.pi * AR_t) / (2 + math.sqrt(4 + AR_t**2))

    NP = 0.25 + (S_tail * l_H) / (S_w * MAC) * (a_t / a_w) * (1 - 2 * a_w / (math.pi * AR_w))
    return (NP - SM) * MAC

def visualize_stl(stl_path):
    if os.path.exists(stl_path):
        mesh = pv.read(stl_path)
        plotter = pv.Plotter(title="Tail Sizing Result")
        plotter.add_mesh(mesh, color="lightblue", show_edges=True)
        plotter.show()

if __name__ == "__main__":
    main()