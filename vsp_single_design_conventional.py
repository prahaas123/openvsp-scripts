import os
import pyvista as pv
import csv
import glob
import subprocess

vsp_exe = r"C:\Program Files\OpenVSP-3.47.0\vsp.exe"

wing_span_res = 20
wing_chord_res = 50
velocity = 10 # m/s
alphas = list(range(-5, 15)) # degrees AoA

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
    "chord": 0.132,
    "l_H": 0.379,        # [m] Tail Moment Arm (Distance from CG to Tail AC)
    "airfoil": "0012",
    "span": 0.45,
    "alpha": -6.85
}

vtail_params = {
    "chord": 0.142,
    "airfoil": "0012",
    "span": 0.18,
    "taper": 0.75,
    "sweep": 20.0
}

def main():
    Sref = wing_params["root_chord"] * wing_params["span"]
    bref = wing_params["span"]
    cref = wing_params["root_chord"]

    # Aero sweep
    stl_path, vsp3_path = generate_wing_and_tail("plane")
    visualize_stl(stl_path)

    CL, CD, Cm = vsp_sweep(vsp3_path, Sref, bref, cref)
    aero_results = zip(alphas, CL, CD, Cm)

    for filename in glob.glob("plane*"):
        try:
            os.remove(filename)
        except OSError:
            pass

    aero_filename = "cfd_sweep.csv"
    aero_headers = ["Alpha_deg", "CL", "CD", "Cm"]
    with open(aero_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(aero_headers)
        writer.writerows(aero_results)

    # Stability sweep
    stl_path, vsp3_path = generate_wing_and_tail("plane")
    vsp_stability(vsp3_path, Sref, bref, cref)

    read_stability("plane.stab")

    for filename in glob.glob("plane*"):
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
        f'    SetParmVal( htail_id, "TotalSpan",      "WingGeom", {htail_params["span"]} );',
        f'    SetParmVal( htail_id, "Root_Chord",     "XSec_1",   {htail_params["chord"]} );',
        f'    SetParmVal( htail_id, "Tip_Chord",      "XSec_1",   {htail_params["chord"]} );',
        '    SetParmVal( htail_id, "Sweep",           "XSec_1",   0.0 );',
        '    SetParmVal( htail_id, "Dihedral",        "XSec_1",   0.0 );',
        '    SetParmVal( htail_id, "Twist",           "XSec_1",   0.0 );',
        f'    SetParmVal( htail_id, "SectTess_U",     "XSec_1",   {wing_span_res}.0 );',
        f'    SetParmVal( htail_id, "Tess_W",         "Shape",    {wing_chord_res}.0 );',
        f'    SetParmVal( htail_id, "X_Rel_Location", "XForm",    {htail_params["l_H"]} );',
        f'    SetParmVal( htail_id, "Y_Rel_Rotation", "XForm",    {htail_params["alpha"]} );',
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
        f'    SetParmVal( vtail_id, "TotalSpan",      "WingGeom", {vtail_params["span"]} );',
        f'    SetParmVal( vtail_id, "Root_Chord",     "XSec_1",   {vtail_params["chord"]} );',
        f'    SetParmVal( vtail_id, "Tip_Chord",      "XSec_1",   {vtail_tip} );',
        f'    SetParmVal( vtail_id, "Sweep",          "XSec_1",   {vtail_params["sweep"]} );',
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

def vsp_sweep(vsp3_path, Sref, bref, cref):
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
        f'    Print( "--- Running Meshing ---" );',
        f'    ExecAnalysis( "VSPAEROComputeGeometry" );',
        f'    SetAnalysisInputDefaults( "VSPAEROSweep" );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Sref",           {darr(Sref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "cref",           {darr(cref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "bref",           {darr(bref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaStart",     {darr(float(alphas[0]))}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaEnd",       {darr(float(alphas[-1]))}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "AlphaNpts",      {iarr(len(alphas))}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "MachStart",      {darr(mach)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "MachNpts",       {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Vinf",           {darr(velocity)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Xcg",            {darr(x_cg)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "WakeNumIter",    {iarr(15)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "NCPU",           {iarr(8)}, 0 );',
        f'    Print( "--- Running Aero Sweep ---" );',
        f'    ExecAnalysis( "VSPAEROSweep" );',
        "}",
    ]

    script_path = "plane_sweep.vspscript"
    with open(script_path, 'w') as f:
        f.write("\n".join(script_lines))

    print(f"--- Running aero sweep ({script_path}) ---")
    subprocess.run([vsp_exe, "-script", script_path], check=True)
    os.remove(script_path)

    CL, CD, Cm = parse_polar("plane.polar")
    return CL, CD, Cm

def vsp_stability(vsp3_path, Sref, bref, cref):
    mach = velocity / 343.0

    script_lines = [
        "void main() {",
        f'    ReadVSPFile( "{vsp3_path}" );',
        f'    SetAnalysisInputDefaults( "VSPAEROComputeGeometry" );',
        f'    array< int > thick_set = GetIntAnalysisInput( "VSPAEROComputeGeometry", "GeomSet" );',
        f'    array< int > thin_set = GetIntAnalysisInput( "VSPAEROComputeGeometry", "ThinGeomSet" );',
        f'    thick_set[0] = ( SET_TYPE::SET_NONE );',
        f'    thin_set[0] = ( SET_TYPE::SET_ALL );',
        f'    AutoGroupVSPAEROControlSurfaces();',
        f'    Update();',
        f'    string container_id = FindContainer( "VSPAEROSettings", 0 );',
        f'    array<string> geoms = FindGeoms();',
        f'    string htail_id = geoms[1];',
        f'    string elev_cs_id = GetSubSurf( htail_id, 0 );',
        f'    string elev_parm = FindParm( container_id, "Surf_" + elev_cs_id + "_1_Gain", "HorizontalTail_SS_CONT_0" );',
        f'    SetParmVal( elev_parm, -1.0 );',
        f'    Update();',
        f'    Print( "--- Running Meshing (stability) ---" );',
        f'    ExecAnalysis( "VSPAEROComputeGeometry" );',
        f'    SetAnalysisInputDefaults( "VSPAEROSweep" );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Sref",           {darr(Sref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "cref",           {darr(cref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "bref",           {darr(bref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaStart",     {darr(0.0)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaEnd",       {darr(0.0)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "AlphaNpts",      {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "MachStart",      {darr(mach)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "MachNpts",       {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Vinf",           {darr(100.0)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Xcg",            {darr(x_cg)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "WakeNumIter",    {iarr(15)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "UnsteadyType",   {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Rho",            {darr(1.225)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "NCPU",           {iarr(8)}, 0 );',
        f'    Print( "--- Running Stability Sweep ---" );',
        f'    ExecAnalysis( "VSPAEROSweep" );',
        "}",
    ]

    script_path = "plane_stab.vspscript"
    with open(script_path, 'w') as f:
        f.write("\n".join(script_lines))

    print(f"--- Running stability sweep ({script_path}) ---")
    subprocess.run([vsp_exe, "-script", script_path], check=True)
    os.remove(script_path)

def parse_polar(polar_path):
    CL, CD, Cm = [], [], []
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
                col_cm = tokens.index('CMytot')
                continue
            if col_cl is None:
                continue
            try:
                CL.append(float(tokens[col_cl]))
                CD.append(float(tokens[col_cd]))
                Cm.append(float(tokens[col_cm]))
            except (ValueError, IndexError):
                continue

    return CL, CD, Cm

def read_stability(stab_path, output_file="vsp_derivatives.csv"): 
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
        'CL_de':   get('CL',  col_elevator),  # CL  wrt elevator
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
 
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Derivative", "Value"])
        for key, val in vsp_dict.items():
            writer.writerow([key, val])
 
    print(f"Successfully generated {output_file} formatted for the 6-DOF simulator.")
    
# AngelScript array literal helpers
def iarr(v):
    return f"array<int> = {{{v}}}"

def darr(v):
    return f"array<double> = {{{v}}}"

if __name__ == "__main__":
    main()