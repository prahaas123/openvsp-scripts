import os
import pyvista as pv
import math
import uuid
import glob
import subprocess
 
vsp_exe = r"C:\Program Files\OpenVSP-3.47.0\vsp.exe"
 
wing_span_res = 20
wing_chord_res = 50
velocity = 10 # m/s
alpha = 0 # degrees AoA
SM = 0.10  # Desired Static Margin
 
airfoil_file = r"C:\Users\kprah\Desktop\Prahaas\WatArrow\CFD Automation\Airfoils\goe322.dat"
 
moment_tolerance = 0.05
tail_sizing_iterations = 10
 
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
    "V_H": 0.10,           # Horizontal Tail Volume Coefficient (Typical: 0.4 - 0.6)
    "l_H": 0.65,          # [m] Tail Moment Arm (Distance from CG to Tail LE)
    "airfoil": "0012",
    "aspect_ratio": 3.0
}
 
vtail_params = {
    "V_V": 0.05,          # Vertical Tail Volume Coefficient (Typical: 0.03 - 0.06)
    "airfoil": "0012",
    "aspect_ratio": 1.5
} 
 
def main():
    Sref = wing_params["span"] * wing_params["root_chord"]
    cref = wing_params["root_chord"]
    x_cg = 0.25 * wing_params["root_chord"]
 
    # Generate just main wing to get moment
    name = f"wing_{uuid.uuid4().hex[:8]}"
    generate_wing(name)
    plain_moment = get_moment(name, x_cg, Sref, cref)
    print(f"Pitching Moment Coefficient (Cmy) at alpha={alpha} deg: {plain_moment:.6f} (without a tail)")
 
    # Generate initial tail geometry
    S_tail = (htail_params["V_H"] * wing_params["span"] * wing_params["span"] * wing_params["root_chord"]) / htail_params["l_H"]
    b_tail = math.sqrt(S_tail * htail_params["aspect_ratio"])
    vtail_height = math.sqrt((vtail_params["V_V"] * wing_params["span"] * wing_params["span"] * wing_params["root_chord"]) / htail_params["l_H"])
 
    tail_name = f"plane_{uuid.uuid4().hex[:8]}"
    new_x_cg = calc_cg(S_tail)
    generate_wing_and_tail(tail_name, b_tail, 0.0, vtail_height)
    moment = get_moment(tail_name, new_x_cg, Sref, cref)
 
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
        generate_wing_and_tail(tail_name, b_tail, i_new, vtail_height)
        m_new = get_moment(tail_name, new_x_cg, Sref, cref)
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
            generate_wing_and_tail(tail_name, b_tail, i_new, vtail_height)
            m_new = get_moment(tail_name, x_cg, Sref, cref)
            print(f"Iter {iteration}: Tail Alpha = {i_new:.2f} deg -> CMy = {m_new:.6f}")
        else:
            print("WARNING: Max iterations reached without full convergence.")
 
    if success:
        tail_chord = b_tail / htail_params["aspect_ratio"]
        print(f"Final Tail Alpha for Trim: {i_new:.4f} degrees with CMy = {m_new:.6f}")
        print(f"Horizontal tail dimensions: Chord: {tail_chord:.3f}m, Span: {b_tail:.3f}m")
        print(f"Final CG location: {new_x_cg:.3f}m from LE, for a Static Margin of {SM*100:.1f}%")
        generate_wing_and_tail("plane_final", b_tail, i_new, vtail_height)
        visualize_stl("plane_final.stl")

def generate_wing(wing_name):
    airfoil_fwd = airfoil_file.replace("\\", "/")
    tip_chord   = wing_params["root_chord"] * wing_params["taper"]
 
    script_lines = [
        "void main() {",
        "    VSPCheckSetup();",
        "    ClearVSPModel();",
        '    string wing_id = AddGeom( "WING" );',
        f'    SetParmVal( wing_id, "TotalSpan",      "WingGeom", {wing_params["span"]} );',
        f'    SetParmVal( wing_id, "Root_Chord",     "XSec_1",   {wing_params["root_chord"]} );',
        f'    SetParmVal( wing_id, "Tip_Chord",      "XSec_1",   {tip_chord} );',
        f'    SetParmVal( wing_id, "Sweep",          "XSec_1",   {wing_params["sweep"]} );',
        f'    SetParmVal( wing_id, "Dihedral",       "XSec_1",   {wing_params["dihedral"]} );',
        f'    SetParmVal( wing_id, "Twist",          "XSec_1",   {wing_params["twist"]} );',
        f'    SetParmVal( wing_id, "Twist_Location", "XSec_1",   1.0 );',
        f'    SetParmVal( wing_id, "SectTess_U",     "XSec_1",   {wing_span_res}.0 );',
        f'    SetParmVal( wing_id, "Tess_W",         "Shape",    {wing_chord_res}.0 );',
        '    string root_surf = GetXSecSurf( wing_id, 0 );',
        '    ChangeXSecShape( root_surf, 0, XS_FILE_AIRFOIL );',
        '    string root_xsec = GetXSec( root_surf, 0 );',
        f'    ReadFileAirfoil( root_xsec, "{airfoil_fwd}" );',
        '    string tip_surf = GetXSecSurf( wing_id, 1 );',
        '    ChangeXSecShape( tip_surf, 1, XS_FILE_AIRFOIL );',
        '    string tip_xsec = GetXSec( tip_surf, 1 );',
        f'    ReadFileAirfoil( tip_xsec, "{airfoil_fwd}" );',
        '    SetSetFlag( wing_id, 1, true );',
        '    Update();',
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
    print(f"VSP file saved: {wing_name}.vsp3")

def generate_wing_and_tail(plane_name, htail_b, htail_alpha, vtail_height):
    airfoil_fwd  = airfoil_file.replace("\\", "/")
    tip_chord    = wing_params["root_chord"] * wing_params["taper"]
    tail_chord  = htail_b / htail_params["aspect_ratio"]

    def naca4(code):
        return int(code[0]) / 100.0, int(code[1]) / 10.0, int(code[2:]) / 100.0
    h_camber, h_cam_loc, h_thick = naca4(htail_params["airfoil"])
    v_camber, v_cam_loc, v_thick = naca4(vtail_params["airfoil"])

    script_lines = [
        "void main() {",
        "    VSPCheckSetup();",
        "    ClearVSPModel();",

        # ---- Main wing ----
        '    string wing_id = AddGeom( "WING" );',
        '    SetGeomName( wing_id, "MainWing" );',
        f'    SetParmVal( wing_id, "TotalSpan",      "WingGeom", {wing_params["span"]} );',
        f'    SetParmVal( wing_id, "Root_Chord",     "XSec_1",   {wing_params["root_chord"]} );',
        f'    SetParmVal( wing_id, "Tip_Chord",      "XSec_1",   {tip_chord} );',
        f'    SetParmVal( wing_id, "Sweep",          "XSec_1",   {wing_params["sweep"]} );',
        f'    SetParmVal( wing_id, "Dihedral",       "XSec_1",   {wing_params["dihedral"]} );',
        f'    SetParmVal( wing_id, "Twist",          "XSec_1",   {wing_params["twist"]} );',
        f'    SetParmVal( wing_id, "Twist_Location", "XSec_1",   1.0 );',
        f'    SetParmVal( wing_id, "Y_Rel_Rotation", "XForm",    {wing_params["alpha"]} );',
        f'    SetParmVal( wing_id, "SectTess_U",     "XSec_1",   {wing_span_res}.0 );',
        f'    SetParmVal( wing_id, "Tess_W",         "Shape",    {wing_chord_res}.0 );',
        '    string w_root_surf = GetXSecSurf( wing_id, 0 );',
        '    ChangeXSecShape( w_root_surf, 0, XS_FILE_AIRFOIL );',
        '    string w_root_xsec = GetXSec( w_root_surf, 0 );',
        f'    ReadFileAirfoil( w_root_xsec, "{airfoil_fwd}" );',
        '    string w_tip_surf = GetXSecSurf( wing_id, 1 );',
        '    ChangeXSecShape( w_tip_surf, 1, XS_FILE_AIRFOIL );',
        '    string w_tip_xsec = GetXSec( w_tip_surf, 1 );',
        f'    ReadFileAirfoil( w_tip_xsec, "{airfoil_fwd}" );',
        '    AddSubSurf( wing_id, SS_CONTROL );',
        '    SetParmVal( wing_id, "SE_Const_Flag",  "SS_Control_1", 1.0 );',
        '    SetParmVal( wing_id, "Length_C_Start", "SS_Control_1", 0.2 );',
        '    SetParmVal( wing_id, "EtaFlag",        "SS_Control_1", 1.0 );',
        '    SetParmVal( wing_id, "EtaStart",       "SS_Control_1", 0.2 );',
        '    SetParmVal( wing_id, "EtaEnd",         "SS_Control_1", 0.8 );',
        '    SetSetFlag( wing_id, 1, true );',

        # ---- Horizontal tail ----
        '    string htail_id = AddGeom( "WING" );',
        '    SetGeomName( htail_id, "HorizontalTail" );',
        f'    SetParmVal( htail_id, "TotalSpan",      "WingGeom", {htail_b} );',
        f'    SetParmVal( htail_id, "Root_Chord",     "XSec_1",   {tail_chord} );',
        f'    SetParmVal( htail_id, "Tip_Chord",      "XSec_1",   {tail_chord} );',
        '    SetParmVal( htail_id, "Sweep",           "XSec_1",   0.0 );',
        '    SetParmVal( htail_id, "Dihedral",        "XSec_1",   0.0 );',
        '    SetParmVal( htail_id, "Twist",           "XSec_1",   0.0 );',
        f'    SetParmVal( htail_id, "SectTess_U",     "XSec_1",   {wing_span_res}.0 );',
        f'    SetParmVal( htail_id, "Tess_W",         "Shape",    {wing_chord_res}.0 );',
        f'    SetParmVal( htail_id, "X_Rel_Location", "XForm",    {htail_params["l_H"]} );',
        f'    SetParmVal( htail_id, "Y_Rel_Rotation", "XForm",    {htail_alpha} );',
        f'    SetParmVal( htail_id, "Camber",         "XSecCurve_0", {h_camber} );',
        f'    SetParmVal( htail_id, "CamberLoc",      "XSecCurve_0", {h_cam_loc} );',
        f'    SetParmVal( htail_id, "ThickChord",     "XSecCurve_0", {h_thick} );',
        f'    SetParmVal( htail_id, "Camber",         "XSecCurve_1", {h_camber} );',
        f'    SetParmVal( htail_id, "CamberLoc",      "XSecCurve_1", {h_cam_loc} );',
        f'    SetParmVal( htail_id, "ThickChord",     "XSecCurve_1", {h_thick} );',
        '    AddSubSurf( htail_id, SS_CONTROL );',
        '    SetParmVal( htail_id, "SE_Const_Flag",  "SS_Control_1", 1.0 );',
        '    SetParmVal( htail_id, "Length_C_Start", "SS_Control_1", 0.25 );',
        '    SetParmVal( htail_id, "EtaFlag",        "SS_Control_1", 1.0 );',
        '    SetParmVal( htail_id, "EtaStart",       "SS_Control_1", 0.0 );',
        '    SetParmVal( htail_id, "EtaEnd",         "SS_Control_1", 0.8 );',
        '    SetSetFlag( htail_id, 1, true );',

        # ---- Vertical tail ----
        '    string vtail_id = AddGeom( "WING" );',
        '    SetGeomName( vtail_id, "VerticalTail" );',
        '    SetParmVal( vtail_id, "Sym_Planar_Flag", "Sym",      0.0 );',
        f'    SetParmVal( vtail_id, "TotalSpan",      "WingGeom", {vtail_height} );',
        f'    SetParmVal( vtail_id, "Root_Chord",     "XSec_1",   {tail_chord} );',
        f'    SetParmVal( vtail_id, "Tip_Chord",      "XSec_1",   {tail_chord} );',
        f'    SetParmVal( vtail_id, "Sweep",          "XSec_1",   {0.0} );',
        f'    SetParmVal( vtail_id, "Dihedral",        "XSec_1",   0.0 );',
        f'    SetParmVal( vtail_id, "Twist",           "XSec_1",   0.0 );',
        f'    SetParmVal( vtail_id, "SectTess_U",     "XSec_1",   {wing_span_res}.0 );',
        f'    SetParmVal( vtail_id, "Tess_W",         "Shape",    {wing_chord_res}.0 );',
        f'    SetParmVal( vtail_id, "X_Rel_Location", "XForm",    {htail_params["l_H"]} );',
        f'    SetParmVal( vtail_id, "X_Rel_Rotation",  "XForm",    90.0 );',
        f'    SetParmVal( vtail_id, "Camber",         "XSecCurve_0", {v_camber} );',
        f'    SetParmVal( vtail_id, "CamberLoc",      "XSecCurve_0", {v_cam_loc} );',
        f'    SetParmVal( vtail_id, "ThickChord",     "XSecCurve_0", {v_thick} );',
        f'    SetParmVal( vtail_id, "Camber",         "XSecCurve_1", {v_camber} );',
        f'    SetParmVal( vtail_id, "CamberLoc",      "XSecCurve_1", {v_cam_loc} );',
        f'    SetParmVal( vtail_id, "ThickChord",     "XSecCurve_1", {v_thick} );',
        f'    AddSubSurf( vtail_id, SS_CONTROL );',
        f'    SetParmVal( vtail_id, "SE_Const_Flag",  "SS_Control_1", 1.0 );',
        f'    SetParmVal( vtail_id, "Length_C_Start", "SS_Control_1", 0.35 );',
        f'    SetParmVal( vtail_id, "EtaFlag",        "SS_Control_1", 1.0 );',
        f'    SetParmVal( vtail_id, "EtaStart",       "SS_Control_1", 0.2 );',
        f'    SetParmVal( vtail_id, "EtaEnd",         "SS_Control_1", 0.8 );',
        f'    SetSetFlag( vtail_id, 1, true );',

        f'    Update();',
        f'    WriteVSPFile( "{plane_name}.vsp3", SET_ALL );',
        f'    ExportFile( "{plane_name}.stl", 0, EXPORT_STL );',
        "}",
    ]

    script_path = f"{plane_name}_geom.vspscript"
    with open(script_path, 'w') as f:
        f.write("\n".join(script_lines))

    print(f"--- Running geometry generation ({script_path}) ---")
    subprocess.run([vsp_exe, "-script", script_path], check=True)
    os.remove(script_path)

    stl_path  = f"{plane_name}.stl"
    vsp3_path = f"{plane_name}.vsp3"
    print(f"STL generated: {stl_path}")
    print(f"VSP file saved: {vsp3_path}")
    return stl_path, vsp3_path

def visualize_stl(stl_path):
    if os.path.exists(stl_path):
        mesh = pv.read(stl_path)
        plotter = pv.Plotter(title="Conventional Plane")
        plotter.add_mesh(mesh, color="lightblue", show_edges=True, smooth_shading=True)
        plotter.add_axes()
        plotter.add_floor(face='-z', i_resolution=10, j_resolution=10, color='gray', opacity=0.2)
        print("Opening PyVista window...")
        plotter.show()
    else:
        print("Error: STL not found.")

def get_moment(model_name, x_cg, Sref, cref):
    vsp3_path = f"{model_name}.vsp3"
    mach = velocity / 343.0
 
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
        f'    ExecAnalysis( "VSPAEROComputeGeometry" );',
        f'    SetAnalysisInputDefaults( "VSPAEROSweep" );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Sref",           {darr(Sref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "cref",           {darr(cref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "bref",           {darr(wing_params["span"])}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Xcg",            {darr(x_cg)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaStart",     {darr(float(alpha))}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaEnd",       {darr(float(alpha))}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "AlphaNpts",      {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "MachStart",      {darr(mach)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "MachNpts",       {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Vinf",           {darr(100.0)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "WakeNumIter",    {iarr(15)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "NCPU",           {iarr(8)}, 0 );',
        f'    SetStringAnalysisInput( "VSPAEROSweep", "RedirectFile", array<string> = {{"{vsp3_path}_log.txt"}}, 0 );',
        f'    ExecAnalysis( "VSPAEROSweep" );',
        "}",
    ]
 
    script_path = f"{model_name}_moment.vspscript"
    with open(script_path, 'w') as f:
        f.write("\n".join(script_lines))
 
    subprocess.run([vsp_exe, "-script", script_path], check=True)
    os.remove(script_path)
 
    # Parse CMy from polar file
    polar_path = f"{model_name}.polar"
    moment = parse_moment(polar_path)
 
    # Clean up all files for this model name
    for filename in glob.glob(f"{model_name}*"):
        try:
            os.remove(filename)
        except OSError:
            pass
 
    return moment

def parse_moment(polar_path):
    col_cm = None
    with open(polar_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            tokens = stripped.split()
            if tokens[0] == 'Beta':
                col_cm = tokens.index('CMytot')
                continue
            if col_cm is None:
                continue
            try:
                return float(tokens[col_cm])
            except (ValueError, IndexError):
                continue
    raise ValueError(f"CMytot not found in {polar_path}")

def calc_cg(S_tail, SM=SM):
    b, c_r, taper = wing_params["span"], wing_params["root_chord"], wing_params["taper"]
    l_H, AR_t = htail_params["l_H"], htail_params["aspect_ratio"]
    AR_w = b**2 / (0.5 * (c_r + taper * c_r) * b)
    sweep_qc = math.atan(math.tan(math.radians(wing_params["sweep"])) - (4/AR_w) * (0.25 * (1 - taper) / (1 + taper)))

    S_w = 0.5 * (c_r + taper * c_r) * b
    MAC = (2/3) * c_r * (1 + taper + taper**2) / (1 + taper)
    AR_w = b**2 / S_w
    sweep_qc = math.atan(math.tan(math.radians(wing_params["sweep"])) - (4/AR_w) * (0.25 * (1 - taper) / (1 + taper)))

    a_w = (2 * math.pi * AR_w) / (2 + math.sqrt(4 + AR_w**2 * (1 + math.tan(sweep_qc)**2)))
    a_t = (2 * math.pi * AR_t) / (2 + math.sqrt(4 + AR_t**2))

    NP = 0.25 + (S_tail * l_H) / (S_w * MAC) * (a_t / a_w) * (1 - 2 * a_w / (math.pi * AR_w))
    x_cg = (NP - SM) * MAC

    print(f"CG location: {x_cg:.4f} m from LE for Static Margin of {SM*100:.1f}%")
    return x_cg

# AngelScript array literal helpers
def iarr(v):
    return f"array<int> = {{{v}}}"
 
def darr(v):
    return f"array<double> = {{{v}}}"
 
if __name__ == "__main__":
    main()