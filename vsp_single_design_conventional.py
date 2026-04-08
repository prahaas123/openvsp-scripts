import os
import pyvista as pv
import csv
import glob
import subprocess
import numpy as np

vsp_exe = r"C:\Program Files\OpenVSP-3.47.0\vsp.exe"

wing_span_res = 20
wing_chord_res = 50
velocities = [10] # m/s
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
    "alpha": 0.9423
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
    stl_path, vsp3_path = generate_wing_and_tail("plane")
    visualize_stl(stl_path)

    # Aero sweep
    csv_exists = False
    for v in velocities:
        print(f"\n=== Running VSP Aero Sweep at {v} m/s ===")
        CL, CD, CDi, Cm = vsp_sweep(vsp3_path, v, Sref, bref, cref)
        lod_filename = "plane.lod"
        cl_data = read_lift_distribution(lod_filename)
        first_aoa_key = list(cl_data.keys())[0]
        span_locations = cl_data[first_aoa_key]['SoverB']
        aero_headers = ["Velocity", "Alpha_deg", "CL", "CD", "Cm"] + [f"Cl_span_{loc:.4f}" for loc in span_locations] + ["Oswald_efficiency"]

        # Combine sweep results with local spanwise results
        aero_results = []
        for i, alpha in enumerate(alphas):
            e = [compute_oswald(CL[i], CDi[i], Sref, bref)]
            base_row = [v, alpha, CL[i], CD[i], Cm[i]]
            # Safely find the corresponding AoA in the parsed dictionary to avoid errors
            closest_aoa_key = min(cl_data.keys(), key=lambda k: abs(k - alpha))
            spanwise_cls = cl_data[closest_aoa_key]['Cl']
            full_row = base_row + spanwise_cls + e
            aero_results.append(full_row)

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
        stl_path, vsp3_path = generate_wing_and_tail("plane")
        vsp_stability(vsp3_path, v, Sref, bref, cref)
        stab_dict = read_stability("plane.stab")
        sm = get_sm("plane.stab", cg=x_cg, mac=cref)
        
        stability_filename = "stability.csv"
        stab_headers = ["Velocity"] + list(stab_dict.keys()) + ["StaticMargin"]
        stab_columns = [v] + list(stab_dict.values()) + [sm]
        
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

def vsp_sweep(vsp3_path, velocity, Sref, bref, cref):
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
        f'    SetStringAnalysisInput( "VSPAEROSweep", "RedirectFile", array<string> = {{"{vsp3_path}_log.txt"}}, 0 );',
        f'    ExecAnalysis( "VSPAEROSweep" );',
        "}",
    ]

    script_path = "plane_sweep.vspscript"
    with open(script_path, 'w') as f:
        f.write("\n".join(script_lines))

    print(f"--- Running aero sweep ({script_path}) ---")
    subprocess.run([vsp_exe, "-script", script_path], check=True)
    os.remove(script_path)

    CL, CD, CDi, Cm = parse_polar("plane.polar")
    return CL, CD, CDi, Cm

def vsp_stability(vsp3_path, velocity, Sref, bref, cref):
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
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Vinf",           {darr(velocity)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Xcg",            {darr(x_cg)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "WakeNumIter",    {iarr(15)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "UnsteadyType",   {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Rho",            {darr(1.225)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "NCPU",           {iarr(8)}, 0 );',
        f'    Print( "--- Running Stability Sweep ---" );',
        f'    SetStringAnalysisInput( "VSPAEROSweep", "RedirectFile", array<string> = {{"{vsp3_path}_log.txt"}}, 0 );',
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
    
def get_sm(stab_path, cg=x_cg, mac=1.0):
    with open(stab_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            
            if stripped.startswith('X_np'):
                parts = stripped.split()
                try:
                    np_val = float(parts[1])
                    sm = (np_val - cg) / mac
                    return sm
                    
                except (ValueError, IndexError):
                    return float('nan')
                    
    return float('nan')
    
# AngelScript array literal helpers
def iarr(v):
    return f"array<int> = {{{v}}}"

def darr(v):
    return f"array<double> = {{{v}}}"

if __name__ == "__main__":
    main()