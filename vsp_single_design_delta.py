import os
import pyvista as pv
import csv
import glob
import subprocess

vsp_exe = r"C:\Program Files\OpenVSP-3.47.0\vsp.exe"

wing_span_res = 20
wing_chord_res = 50
velocity = 20.0 # m/s (for stability sweep)
speeds = list(range(5, 45, 5)) # m/s
alphas = list(range(-5, 15)) # degrees AoA

airfoil_file = r"C:\Users\kprah\Desktop\Prahaas\WatArrow\CFD Automation\Airfoils\dae21.dat"

root_chord = 0.2556
taper_ratio = 0.3377
sweep = 42.2426
dihedral = 5.0
twist = -1.8498
span = 0.9034
x_cg = 0.2082
elevon_length = 0.1         # % of chord
elevon_start = 0.2          # % of wingspan
elevon_end = 0.8            # % of wingspan

def main():
    # Aero sweep
    stl_path, vsp3_path = generate_wing("wing")
    visualize_stl(stl_path)

    Mach, AoA, CL, CD, Cm = vsp_sweep(vsp3_path)
    speeds = [round(m * 343.0, 1) for m in Mach]
    aero_results = zip(speeds, AoA, CL, CD, Cm)

    for filename in glob.glob("wing*"):
        try:
            os.remove(filename)
        except OSError:
            pass

    aero_filename = "cfd_sweep.csv"
    aero_headers = ["Speed_m/s", "Alpha_deg", "CL", "CD", "Cm"]
    with open(aero_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(aero_headers)
        writer.writerows(aero_results)

    # # Stability sweep
    stl_path, vsp3_path = generate_wing("wing")
    vsp_stability(vsp3_path)

    read_stability("wing.stab")

    for filename in glob.glob("wing*"):
        try:
            os.remove(filename)
        except OSError:
            pass

def generate_wing(wing_name):
    tip_chord = root_chord * taper_ratio
    airfoil_fwd = airfoil_file.replace("\\", "/")

    script_lines = [
        "void main() {",
        "    VSPCheckSetup();",
        "    ClearVSPModel();",
        f'    string wing_id = AddGeom( "WING" );',
        f'    SetParmVal( wing_id, "TotalSpan",      "WingGeom", {span} );',
        f'    SetParmVal( wing_id, "Root_Chord",     "XSec_1",   {root_chord} );',
        f'    SetParmVal( wing_id, "Tip_Chord",      "XSec_1",   {tip_chord} );',
        f'    SetParmVal( wing_id, "Sweep",          "XSec_1",   {sweep} );',
        f'    SetParmVal( wing_id, "Dihedral",       "XSec_1",   {dihedral} );',
        f'    SetParmVal( wing_id, "Twist",          "XSec_1",   {twist} );',
        f'    SetParmVal( wing_id, "Twist_Location", "XSec_1",   1.0 );',
        f'    SetParmVal( wing_id, "SectTess_U",     "XSec_1",   {wing_span_res}.0 );',
        f'    SetParmVal( wing_id, "Tess_W",         "Shape",    {wing_chord_res}.0 );',
        f'    string root_surf = GetXSecSurf( wing_id, 0 );',
        f'    ChangeXSecShape( root_surf, 0, XS_FILE_AIRFOIL );',
        f'    string root_xsec = GetXSec( root_surf, 0 );',
        f'    ReadFileAirfoil( root_xsec, "{airfoil_fwd}" );',
        f'    string tip_surf = GetXSecSurf( wing_id, 1 );',
        f'    ChangeXSecShape( tip_surf, 1, XS_FILE_AIRFOIL );',
        f'    string tip_xsec = GetXSec( tip_surf, 1 );',
        f'    ReadFileAirfoil( tip_xsec, "{airfoil_fwd}" );',
        f'    SetSetFlag( wing_id, 1, true );',
        f'    AddSubSurf( wing_id, SS_CONTROL );',
        f'    SetParmVal( wing_id, "SE_Const_Flag",  "SS_Control_1", 1.0 );',
        f'    SetParmVal( wing_id, "Length_C_Start", "SS_Control_1", {elevon_length} );',
        f'    SetParmVal( wing_id, "EtaFlag",        "SS_Control_1", 1.0 );',
        f'    SetParmVal( wing_id, "EtaStart",       "SS_Control_1", {elevon_start} );',
        f'    SetParmVal( wing_id, "EtaEnd",         "SS_Control_1", {elevon_end} );',
        f'    Update();',
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

    stl_path = f"{wing_name}.stl"
    vsp3_path = f"{wing_name}.vsp3"
    print(f"STL generated: {stl_path}")
    print(f"VSP file saved: {vsp3_path}")
    return stl_path, vsp3_path

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

def vsp_sweep(vsp3_path):
    Sref = 0.5 * (root_chord + root_chord * taper_ratio) * span
    bref = span
    cref = root_chord
    machs = [v / 343.0 for v in speeds]

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
        f'    Print( "--- Running Meshing ---" );',
        f'    ExecAnalysis( "VSPAEROComputeGeometry" );',
        f'    SetAnalysisInputDefaults( "VSPAEROSweep" );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Sref",        {darr(Sref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "cref",        {darr(cref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "bref",        {darr(bref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaStart",  {darr(float(alphas[0]))}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaEnd",    {darr(float(alphas[-1]))}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "AlphaNpts",   {iarr(len(alphas))}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "MachStart",   {darr(float(machs[0]))}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "MachEnd",     {darr(float(machs[-1]))}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "MachNpts",    {iarr(len(machs))}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Vinf",        {darr(20.0)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Xcg",         {darr(x_cg)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "WakeNumIter", {iarr(15)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "NCPU",        {iarr(8)}, 0 );',
        f'    Print( "--- Running Aero Sweep ---" );',
        f'    SetStringAnalysisInput( "VSPAEROSweep", "RedirectFile", array<string> = {{"{vsp3_path}_log.txt"}}, 0 );',
        f'    ExecAnalysis( "VSPAEROSweep" );',
        "}",
    ]

    script_path = "wing_sweep.vspscript"
    with open(script_path, 'w') as f:
        f.write("\n".join(script_lines))

    print(f"--- Running aero sweep ({script_path}) ---")
    subprocess.run([vsp_exe, "-script", script_path], check=True)
    os.remove(script_path)

    Mach, AoA, CL, CD, Cm = parse_polar("wing.polar")
    return Mach, AoA, CL, CD, Cm

def vsp_stability(vsp3_path):
    Sref = 0.5 * (root_chord + root_chord * taper_ratio) * span
    bref = span
    cref = root_chord
    mach = velocity / 343.0

    script_lines = [
        "void main() {",
        f'    ReadVSPFile( "{vsp3_path}" );',
        f'    SetAnalysisInputDefaults( "VSPAEROComputeGeometry" );',
        f'    array< int > thick_set = GetIntAnalysisInput( "VSPAEROComputeGeometry", "GeomSet" );',
        f'    array< int > thin_set = GetIntAnalysisInput( "VSPAEROComputeGeometry", "ThinGeomSet" );',
        f'    thick_set[0] = ( SET_TYPE::SET_NONE );',
        f'    thin_set[0] = ( SET_TYPE::SET_ALL );',
        f'    int cs_pitch_id = CreateVSPAEROControlSurfaceGroup();',
        f'    SetVSPAEROControlGroupName( "Pitch", cs_pitch_id );',
        f'    AddAllToVSPAEROControlSurfaceGroup( cs_pitch_id );',
        f'    int cs_roll_id = CreateVSPAEROControlSurfaceGroup();',
        f'    SetVSPAEROControlGroupName( "Roll", cs_roll_id );',
        f'    AddAllToVSPAEROControlSurfaceGroup( cs_roll_id );',
        f'    Update();',
        f'    string container_id = FindContainer( "VSPAEROSettings", 0 );',
        f'    array<string> geoms = FindGeoms();',
        f'    string wid = geoms[0];',
        f'    string cs_id = GetSubSurf( wid, 0 );',
        f'    string group_str = "ControlSurfaceGroup_" + formatInt( cs_pitch_id + 1, "" );',
        f'    string parm_id = FindParm( container_id, "Surf_" + cs_id + "_1_Gain", group_str );',
        f'    SetParmVal( parm_id, -1.0 );',
        f'    Print( "--- Running Meshing (stability) ---" );',
        f'    ExecAnalysis( "VSPAEROComputeGeometry" );',
        f'    SetAnalysisInputDefaults( "VSPAEROSweep" );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Sref",         {darr(Sref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "cref",         {darr(cref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "bref",         {darr(bref)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaStart",   {darr(0.0)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "AlphaEnd",     {darr(0.0)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "AlphaNpts",    {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "MachStart",    {darr(mach)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "MachNpts",     {iarr(1)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Vinf",         {darr(100.0)}, 0 );',
        f'    SetDoubleAnalysisInput( "VSPAEROSweep", "Xcg",          {darr(x_cg)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "WakeNumIter",  {iarr(15)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "UnsteadyType", {iarr(1)}, 0 );',
        f'    SetIntAnalysisInput(    "VSPAEROSweep", "NCPU",         {iarr(8)}, 0 );',
        f'    Print( "--- Running Stability Sweep ---" );',
        f'    SetStringAnalysisInput( "VSPAEROSweep", "RedirectFile", array<string> = {{"{vsp3_path}_log.txt"}}, 0 );',
        f'    ExecAnalysis( "VSPAEROSweep" );',
        "}",
    ]

    script_path = "wing_stab.vspscript"
    with open(script_path, 'w') as f:
        f.write("\n".join(script_lines))

    print(f"--- Running stability sweep ({script_path}) ---")
    subprocess.run([vsp_exe, "-script", script_path], check=True)
    os.remove(script_path)

def parse_polar(polar_path):
    Mach, AoA, CL, CD, Cm = [], [], [], [], []
    col_mach = col_aoa = col_cl = col_cd = col_cm = None
 
    with open(polar_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            tokens = stripped.split()
            if tokens[0] == 'Beta':
                col_mach = tokens.index('Mach')
                col_aoa  = tokens.index('AoA')
                col_cl   = tokens.index('CLtot')
                col_cd   = tokens.index('CDtot')
                col_cm   = tokens.index('CMytot')
                continue
            if col_cl is None:
                continue
            try:
                Mach.append(float(tokens[col_mach]))
                AoA.append(float(tokens[col_aoa]))
                CL.append(float(tokens[col_cl]))
                CD.append(float(tokens[col_cd]))
                Cm.append(float(tokens[col_cm]))
            except (ValueError, IndexError):
                continue
 
    return Mach, AoA, CL, CD, Cm

def read_stability(stab_path, output_file="vsp_derivatives.csv"):
    col_alpha = col_beta = col_p = col_q = col_r = col_pitch = col_roll = None
    rows = {}
 
    with open(stab_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            tokens = stripped.split()
            if tokens[0] == 'Coef':
                col_alpha = tokens.index('Alpha')
                col_beta  = tokens.index('Beta')
                col_p     = tokens.index('p')
                col_q     = tokens.index('q')
                col_r     = tokens.index('r')
                col_pitch = tokens.index('ConGrp_1')
                col_roll  = tokens.index('ConGrp_2')
                continue
            if col_alpha is None:
                continue
            coef = tokens[0]
            if coef in ('CL', 'CS', 'CMl', 'CMm', 'CMn', 'CD'):
                try:
                    rows[coef] = [float(tokens[col_alpha]),
                                  float(tokens[col_beta]),
                                  float(tokens[col_p]),
                                  float(tokens[col_q]),
                                  float(tokens[col_r]),
                                  float(tokens[col_pitch]),
                                  float(tokens[col_roll])]
                except (ValueError, IndexError):
                    pass
 
    vsp_dict = {}
    if 'CL'  in rows: vsp_dict['CL_de']   = rows['CL'] [5]  # CL   wrt ConGrp_1 (pitch)
    if 'CS'  in rows: vsp_dict['CY_beta'] = rows['CS'] [1]  # CS   wrt Beta
    if 'CS'  in rows: vsp_dict['CY_p']    = rows['CS'] [2]  # CS   wrt p
    if 'CS'  in rows: vsp_dict['CY_r']    = rows['CS'] [4]  # CS   wrt r
    if 'CMl' in rows: vsp_dict['Cl_beta'] = rows['CMl'][1]  # CMl  wrt Beta
    if 'CMl' in rows: vsp_dict['Cl_p']    = rows['CMl'][2]  # CMl  wrt p
    if 'CMl' in rows: vsp_dict['Cl_r']    = rows['CMl'][4]  # CMl  wrt r
    if 'CMl' in rows: vsp_dict['Cl_da']   = rows['CMl'][6]  # CMl  wrt ConGrp_2 (roll)
    if 'CMm' in rows: vsp_dict['Cm_q']    = rows['CMm'][3]  # CMm  wrt q
    if 'CMm' in rows: vsp_dict['Cm_de']   = rows['CMm'][5]  # CMm  wrt ConGrp_1 (pitch)
    if 'CMn' in rows: vsp_dict['Cn_beta'] = rows['CMn'][1]  # CMn  wrt Beta
    if 'CMn' in rows: vsp_dict['Cn_p']    = rows['CMn'][2]  # CMn  wrt p
    if 'CMn' in rows: vsp_dict['Cn_r']    = rows['CMn'][4]  # CMn  wrt r
    if 'CMn' in rows: vsp_dict['Cn_da']   = rows['CMn'][6]  # CMn  wrt ConGrp_2 (roll)
 
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Derivative", "Value"])
        for key, val in vsp_dict.items():
            writer.writerow([key, val])
 
    print("Successfully generated " + output_file + " formatted for the 6-DOF simulator.")

# AngelScript array helpers
def iarr(v):
    return f"array<int> = {{{v}}}"

def darr(v):
    return f"array<double> = {{{v}}}"

if __name__ == "__main__":
    main()