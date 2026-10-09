import os
import csv
import itertools
import openvsp as vsp # type: ignore
import math
import uuid
import glob

wing_span_res = 20
wing_chord_res = 50
velocity = 10 # m/s
alpha = 0 # degrees AoA
SM = 0.10  # desired Static Margin

airfoil_file = r"Airfoils\goe322.dat"
output_csv = "tail_sizing_options.csv"

moment_tolerance = 0.05
tail_sizing_iterations = 10

# Every combination of wing config x tail moment arm x tail AR is sized
wing_configs = [
    {
        "name": "baseline",
        "span": 1.02,
        "root_chord": 0.27,
        "taper": 1.0,
        "sweep": 0.0,
        "dihedral": 0.0,
        "twist": 0.0,
        "alpha": 3.0
    },
    {
        "name": "long_span",
        "span": 1.20,
        "root_chord": 0.25,
        "taper": 1.0,
        "sweep": 0.0,
        "dihedral": 0.0,
        "twist": 0.0,
        "alpha": 3.0
    },
]

tail_moment_arms = [0.45, 0.52, 0.60]  # l_H [m]
htail_aspect_ratios = [2.5, 3.0, 3.5]

htail_params = {
    "V_H": 0.6,
    "airfoil": "0012"
}

vtail_params = {
    "V_V": 0.05,
    "airfoil": "0012"
}

csv_fields = [
    "wing_name", "wing_span", "wing_root_chord", "wing_taper", "wing_sweep",
    "wing_dihedral", "wing_twist", "wing_incidence", "wing_area",
    "l_H", "htail_AR", "htail_span", "htail_chord", "vtail_height",
    "tail_incidence", "converged"
]

def main():
    vsp.VSPCheckSetup()

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields)
        writer.writeheader()

        for wing in wing_configs:
            S_w, MAC = wing_geometry(wing)

            # Generate just main wing to get initial moment
            name = f"wing_{uuid.uuid4().hex[:8]}"
            vsp3_path = generate_wing(name, wing)
            plain_moment = get_moment(vsp3_path, 0.25 * MAC, S_w, MAC, wing["span"])
            print(f"[{wing['name']}] Cmy at alpha={alpha} deg: {plain_moment:.6f} (without a tail)")

            for l_H, AR_t in itertools.product(tail_moment_arms, htail_aspect_ratios):
                htail = dict(htail_params, l_H=l_H, aspect_ratio=AR_t)
                print(f"\n=== {wing['name']}: l_H = {l_H:.3f} m, tail AR = {AR_t:.2f} ===")
                result = size_tail(wing, htail)

                writer.writerow({
                    "wing_name": wing["name"],
                    "wing_span": wing["span"],
                    "wing_root_chord": wing["root_chord"],
                    "wing_taper": wing["taper"],
                    "wing_sweep": wing["sweep"],
                    "wing_dihedral": wing["dihedral"],
                    "wing_twist": wing["twist"],
                    "wing_incidence": wing["alpha"],
                    "wing_area": round(S_w, 5),
                    "l_H": l_H,
                    "htail_AR": AR_t,
                    **result
                })
                f.flush()  # keep finished rows if a later run crashes

    print(f"\nResults written to {output_csv}")

def size_tail(wing, htail):
    S_w, MAC = wing_geometry(wing)

    # Generate initial tail geometry
    S_tail = htail["V_H"] * S_w * MAC / htail["l_H"]
    b_tail = math.sqrt(S_tail * htail["aspect_ratio"])
    tail_chord = b_tail / htail["aspect_ratio"]

    # Vertical tail shares the horizontal tail chord, so height = S_V / chord
    S_vtail = vtail_params["V_V"] * S_w * wing["span"] / htail["l_H"]
    vtail_height = S_vtail / tail_chord
    x_cg = calc_cg(S_tail, wing, htail)

    def trim_moment(tail_alpha):
        tail_name = f"plane_{uuid.uuid4().hex[:8]}"
        vsp3_path = generate_wing_and_tail(tail_name, wing, htail, b_tail, tail_alpha, vtail_height)
        return get_moment(vsp3_path, x_cg, S_w, MAC, wing["span"])

    # Tail incidence sizing loop (secant method)
    i_old = 0.0
    m_old = trim_moment(i_old)
    success = False
    print(f"Iter 0: Tail Alpha = {i_old:.2f} deg -> CMy = {m_old:.6f}")

    if abs(m_old) < moment_tolerance:
        print("Aircraft is naturally trimmed!")
        i_new = i_old
        m_new = m_old
        success = True
    else:
        i_new = -2.0  # Initial guess
        m_new = trim_moment(i_new)
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

            m_new = trim_moment(i_new)
            print(f"Iter {iteration}: Tail Alpha = {i_new:.2f} deg -> CMy = {m_new:.6f}")
        else:
            if abs(m_new) < moment_tolerance:
                success = True
            else:
                print("WARNING: Max iterations reached without full convergence.")

    print(f"Tail Alpha: {i_new:.4f} deg, CMy = {m_new:.6f}, "
          f"HTail chord {tail_chord:.3f} m x span {b_tail:.3f} m, CG {x_cg:.3f} m from LE")

    return {
        "htail_span": round(b_tail, 5),
        "htail_chord": round(tail_chord, 5),
        "vtail_height": round(vtail_height, 5),
        "tail_incidence": round(i_new, 4),
        "converged": success
    }

def wing_geometry(wing):
    c_r, taper = wing["root_chord"], wing["taper"]
    S_w = 0.5 * (c_r + taper * c_r) * wing["span"]
    MAC = (2/3) * c_r * (1 + taper + taper**2) / (1 + taper)
    return S_w, MAC

def add_main_wing(wing):
    wid = vsp.AddGeom("WING", "")
    tip_chord = wing["root_chord"] * wing["taper"]

    vsp.SetParmVal(wid, "TotalSpan", "WingGeom", wing["span"])
    vsp.SetParmVal(wid, "Root_Chord", "XSec_1", wing["root_chord"])
    vsp.SetParmVal(wid, "Tip_Chord", "XSec_1", tip_chord)
    vsp.SetParmVal(wid, "Sweep", "XSec_1", wing["sweep"])
    vsp.SetParmVal(wid, "Dihedral", "XSec_1", wing["dihedral"])
    vsp.SetParmVal(wid, "Twist", "XSec_1", wing["twist"])
    vsp.SetParmVal(wid, "Twist_Location", "XSec_1", 1.0)
    vsp.SetParmVal(wid, "SectTess_U", "XSec_1", float(wing_span_res))
    vsp.SetParmVal(wid, "Tess_W", "Shape", float(wing_chord_res))

    surf = vsp.GetXSecSurf(wid, 0)
    for i in [0, 1]:
        vsp.ChangeXSecShape(surf, i, vsp.XS_FILE_AIRFOIL)
        vsp.ReadFileAirfoil(vsp.GetXSec(surf, i), airfoil_file)
    vsp.SetSetFlag(wid, 1, True)
    return wid

def generate_wing(wing_name, wing):
    vsp.ClearVSPModel()
    add_main_wing(wing)

    vsp.Update()
    vsp3_path = f"{wing_name}.vsp3"
    vsp.WriteVSPFile(vsp3_path)
    return vsp3_path

def generate_wing_and_tail(plane_name, wing, htail, htail_b, htail_alpha, vtail_height):
    vsp.ClearVSPModel()
    tail_chord = htail_b / htail["aspect_ratio"]

    def naca4(code):
        return int(code[0]) / 100.0, int(code[1]) / 10.0, int(code[2:]) / 100.0

    h_camber, h_cam_loc, h_thick = naca4(htail["airfoil"])
    v_camber, v_cam_loc, v_thick = naca4(vtail_params["airfoil"])

    # Main Wing
    wid = add_main_wing(wing)
    vsp.SetGeomName(wid, "MainWing")
    vsp.SetParmVal(wid, "Y_Rel_Rotation", "XForm", wing["alpha"])

    # Horizontal Tail
    hid = vsp.AddGeom("WING", "")
    vsp.SetGeomName(hid, "HorizontalTail")
    vsp.SetParmVal(hid, "TotalSpan", "WingGeom", htail_b)
    vsp.SetParmVal(hid, "Root_Chord", "XSec_1", tail_chord)
    vsp.SetParmVal(hid, "Tip_Chord", "XSec_1", tail_chord)
    vsp.SetParmVal(hid, "Sweep", "XSec_1", 0.0)
    vsp.SetParmVal(hid, "X_Rel_Location", "XForm", htail["l_H"])
    vsp.SetParmVal(hid, "Y_Rel_Rotation", "XForm", htail_alpha)
    vsp.SetParmVal(hid, "Camber", "XSecCurve_0", h_camber)
    vsp.SetParmVal(hid, "CamberLoc", "XSecCurve_0", h_cam_loc)
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
    vsp.SetParmVal(vid, "X_Rel_Location", "XForm", htail["l_H"])
    vsp.SetParmVal(vid, "X_Rel_Rotation", "XForm", 90.0)
    vsp.SetParmVal(vid, "Camber", "XSecCurve_0", v_camber)
    vsp.SetParmVal(vid, "CamberLoc", "XSecCurve_0", v_cam_loc)
    vsp.SetParmVal(vid, "ThickChord", "XSecCurve_0", v_thick)
    vsp.SetSetFlag(vid, 1, True)

    vsp.Update()
    vsp3_path = f"{plane_name}.vsp3"
    vsp.WriteVSPFile(vsp3_path)
    return vsp3_path

def get_moment(vsp3_path, x_cg, Sref, cref, bref):
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
    vsp.SetDoubleAnalysisInput("VSPAEROSweep", "bref", [bref])
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

def calc_cg(S_tail, wing, htail, SM=SM):
    S_w, MAC = wing_geometry(wing)
    l_H, AR_t = htail["l_H"], htail["aspect_ratio"]
    AR_w = wing["span"]**2 / S_w

    # Simple estimate for lift curves
    a_w = (2 * math.pi * AR_w) / (2 + math.sqrt(4 + AR_w**2))
    a_t = (2 * math.pi * AR_t) / (2 + math.sqrt(4 + AR_t**2))

    NP = 0.25 + (S_tail * l_H) / (S_w * MAC) * (a_t / a_w) * (1 - 2 * a_w / (math.pi * AR_w))
    return (NP - SM) * MAC

if __name__ == "__main__":
    main()
