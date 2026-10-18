import os
import csv
import glob
import math
import uuid
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import qmc

# INPUTS

N_SAMPLES = 256
KEEP_PCTL = (2.0, 98.0)  # percentiles of feasible samples
VELOCITY = 15.0          # [m/s] cruise airspeed

VAR_NAMES = ["root_chord", "taper", "sweep", "washout", "span", "alpha"]
BOUNDS = [
    (0.15, 0.40),    # root_chord [m]
    (0.25, 1.00),    # taper ratio
    (0.0, 45.0),     # LE sweep [deg]
    (-5.0, 0.0),     # washout [deg]
    (0.50, 0.60),    # full span [m]
    (2.0, 8.0),      # alpha [deg]
]

DIHEDRAL = 2.5           # [deg] fixed dihedral value
AIRFOIL_FILE = os.path.join("Airfoils", "mh61.dat")

# Weight Estimate
M_FIXED = 0.265          # [kg]      avionics, battery, propulsion, payload
M_SPAR_PER_M = 0.044     # [kg/m]    spar mass per unit span
M_SKIN_PER_M2 = 0.950    # [kg/m^2]  3D-printed skin + ribs per unit planform area

# Constraint Targets
STATIC_MARGIN = 0.05
CM_TOL = 0.02
AR_MIN = 2.0
AR_MAX = 6.0
TIP_CHORD_MIN = 0.05

# OTHER INPUTS (NO CHANGES REQUIRED)

RHO = 1.225
Q_INF = 0.5 * RHO * VELOCITY ** 2

ALPHA_STEP_DEG = 1.0

SAMPLE_CSV = "dse_samples.csv"
CONSTRAINTS = ["lift", "trim", "ar_min", "ar_max", "tip_chord"]
DERIVED = ["AR", "S_ref", "MAC", "tip_chord", "CL", "CD", "LD", "drag", "lift",
           "weight", "CM_cg", "CM_alpha_cg", "x_np", "x_cg", "feasible"]
FIELDS = VAR_NAMES + DERIVED + [f"v_{c}" for c in CONSTRAINTS]

LO = np.array([b[0] for b in BOUNDS])
HI = np.array([b[1] for b in BOUNDS])
ND = len(VAR_NAMES)

def total_weight(span, s_ref):
    return 9.81 * (M_FIXED + M_SPAR_PER_M * span + M_SKIN_PER_M2 * s_ref)

def planform(root_chord, taper, span):
    tip_chord = root_chord * taper
    s_ref = 0.5 * (root_chord + tip_chord) * span
    ar = span ** 2 / max(s_ref, 1e-12)
    mac = (2.0 / 3.0) * root_chord * (1.0 + taper + taper ** 2) / (1.0 + taper)
    y_mac = (span / 2.0) * (1.0 + 2.0 * taper) / (3.0 * (1.0 + taper))
    return dict(tip_chord=tip_chord, s_ref=s_ref, ar=ar, mac=mac, y_mac=y_mac)

import vsp_optimization as vo  # noqa: E402

def evaluate(x):
    root_chord, taper, sweep, washout, span, alpha = [float(v) for v in x]
    run_id = f"dse_{uuid.uuid4().hex[:8]}"
    try:
        g = planform(root_chord, taper, span)
        s_ref, ar, mac = g["s_ref"], g["ar"], g["mac"]
        weight = total_weight(span, s_ref)
        cd0 = vo.estimate_cd0(mac, s_ref, sweep)

        vsp3_path, _ = vo.generate_wing(run_id, span, root_chord, taper,
                                        sweep, DIHEDRAL, washout)
        r = vo._run_vspaero(vsp3_path, s_ref, span, mac, x_cg=0.0, cd0=cd0,
                            alpha_start=alpha - ALPHA_STEP_DEG,
                            alpha_end=alpha + ALPHA_STEP_DEG, alpha_npts=3)

        a = np.asarray(r["Alpha"], dtype=float)
        cl_m, cl_0, cl_p = [vo._pick(r["CL"], a, alpha + d * ALPHA_STEP_DEG) for d in (-1, 0, 1)]
        cm_m, cm_0, cm_p = [vo._pick(r["CMy"], a, alpha + d * ALPHA_STEP_DEG) for d in (-1, 0, 1)]
        cd_0 = vo._pick(r["CD"], a, alpha)

        d_alpha = 2.0 * math.radians(ALPHA_STEP_DEG)
        cl_alpha, cm_alpha = (cl_p - cl_m) / d_alpha, (cm_p - cm_m) / d_alpha
        if not np.isfinite(cl_alpha) or abs(cl_alpha) < 1e-6:
            raise RuntimeError("degenerate CL_alpha")

        x_np = -mac * (cm_alpha / cl_alpha)
        x_cg = x_np - STATIC_MARGIN * mac
        cm_cg = cm_0 + cl_0 * x_cg / mac
        cm_alpha_cg = cm_alpha + cl_alpha * x_cg / mac

        lift = Q_INF * s_ref * cl_0
        drag = Q_INF * s_ref * cd_0
        ld = cl_0 / cd_0 if cd_0 > 1e-9 else float("nan")

        v = {
            "lift": (weight - lift) / weight,
            "trim": (abs(cm_cg) - CM_TOL) / CM_TOL,
            "ar_min": (AR_MIN - ar) / AR_MIN,
            "ar_max": (ar - AR_MAX) / AR_MAX,
            "tip_chord": (TIP_CHORD_MIN - g["tip_chord"]) / TIP_CHORD_MIN,
        }
        if not all(np.isfinite(list(v.values()))):
            raise RuntimeError("non-finite constraint margin")

        row = dict(zip(VAR_NAMES, [float(t) for t in x]))
        row.update(AR=ar, S_ref=s_ref, MAC=mac, tip_chord=g["tip_chord"],
                   CL=cl_0, CD=cd_0, LD=ld, drag=drag, lift=lift, weight=weight,
                   CM_cg=cm_cg, CM_alpha_cg=cm_alpha_cg, x_np=x_np, x_cg=x_cg,
                   feasible=float(all(val <= 0.0 for val in v.values())))
        row.update({f"v_{k}": val for k, val in v.items()})
        return row

    except Exception as e:
        print(f"  [{run_id}] FAILED: {e}")
        return None
    finally:
        for filename in glob.glob(f"{run_id}*"):
            try:
                os.remove(filename)
            except OSError:
                pass

def append(path, row):
    new = not os.path.isfile(path)
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow({k: round(v, 8) for k, v in row.items()})

def load(path):
    if not os.path.isfile(path):
        return None
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None
    out = {}
    for k in FIELDS:
        out[k] = np.array([float(r.get(k, "nan") or "nan") for r in rows])
    return out

def count_rows(path):
    d = load(path)
    return 0 if d is None else len(d["drag"])

def sample():
    done = count_rows(SAMPLE_CSV)
    if done >= N_SAMPLES:
        print(f"[sample] {done} samples already on disk; skipping.")
        return
    pts = qmc.Sobol(d=ND, scramble=True, seed=1).random(N_SAMPLES)
    pts = LO + pts * (HI - LO)
    print(f"[sample] resuming at {done}/{N_SAMPLES}")
    with vo.make_pool() as pool:
        for i, row in zip(range(done, N_SAMPLES), pool.imap(evaluate, pts[done:])):
            if row is None:
                continue
            append(SAMPLE_CSV, row)
            tag = "FEASIBLE" if row["feasible"] else "infeasible"
            print(f"  [{i + 1}/{N_SAMPLES}] {tag:10s} D={row['drag']:.4f} N  "
                  f"L/D={row['LD']:.2f}  AR={row['AR']:.2f}")

def tighten(d):
    m = d["feasible"] > 0.5
    out = {}
    for j, name in enumerate(VAR_NAMES):
        if m.sum() < 10:
            out[name] = [float(LO[j]), float(HI[j])]
            continue
        lo, hi = np.percentile(d[name][m], KEEP_PCTL)
        out[name] = [float(max(LO[j], lo)), float(min(HI[j], hi))]
    return out

def in_box(d, box):
    m = np.ones(len(d["drag"]), dtype=bool)
    for name in VAR_NAMES:
        m &= (d[name] >= box[name][0]) & (d[name] <= box[name][1])
    return m

def plot_pairs(d, box):
    m = d["feasible"] > 0.5
    fig, ax = plt.subplots(ND - 1, ND - 1, figsize=(2.2 * (ND - 1), 2.2 * (ND - 1)))
    norm = plt.Normalize(*np.percentile(d["drag"][m], [5, 95])) if m.sum() else None
    for r in range(1, ND):
        for c in range(ND - 1):
            a = ax[r - 1, c]
            if c >= r:
                a.axis("off")
                continue
            xn, yn = VAR_NAMES[c], VAR_NAMES[r]
            a.scatter(d[xn][~m], d[yn][~m], s=4, c="0.85", lw=0)
            sc = a.scatter(d[xn][m], d[yn][m], s=8, c=d["drag"][m],
                           cmap="viridis", norm=norm, lw=0)
            a.add_patch(Rectangle((box[xn][0], box[yn][0]),
                                  box[xn][1] - box[xn][0], box[yn][1] - box[yn][0],
                                  fill=False, ec="crimson", lw=1.2))
            a.set_xlim(LO[c], HI[c])
            a.set_ylim(LO[r], HI[r])
            a.tick_params(labelsize=6)
            if r == ND - 1:
                a.set_xlabel(xn, fontsize=8)
            else:
                a.set_xticklabels([])
            if c == 0:
                a.set_ylabel(yn, fontsize=8)
            else:
                a.set_yticklabels([])
    if m.sum():
        fig.colorbar(sc, ax=ax, shrink=0.5, label="drag [N]")
    fig.suptitle("Feasible region and tightened bounds (red)", fontsize=11)

def summarize(d, box):
    m = d["feasible"] > 0.5
    kept = in_box(d, box) & m

    print("\n--- FEASIBILITY ---")
    print(f"samples            : {len(d['drag'])}")
    print(f"feasible fraction  : {m.mean() * 100:.1f}%")
    print(f"after tightening   : {kept.sum() / max(in_box(d, box).sum(), 1) * 100:.1f}% "
          f"({kept.sum() / max(m.sum(), 1) * 100:.0f}% of feasible designs retained)")

    print("\n--- TIGHTENED BOUNDS ---")
    for i, n in enumerate(VAR_NAMES):
        print(f"{n:12s}: [{LO[i]:7.3f}, {HI[i]:7.3f}] -> [{box[n][0]:7.3f}, {box[n][1]:7.3f}]")
    print()

def main():
    sample()
    d = load(SAMPLE_CSV)
    if d is None:
        print("No samples were logged.")
        return

    box = tighten(d)
    summarize(d, box)

    plot_pairs(d, box)
    plt.show()

if __name__ == "__main__":
    main()