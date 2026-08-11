import os
import csv
import glob
import math
import uuid
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import qmc
from scipy.spatial import cKDTree
from scipy.interpolate import RBFInterpolator

from vsp_optimization import (
    VAR_NAMES, BOUNDS, Q_INF, ALPHA_STEP_DEG, STATIC_MARGIN,
    CM_TOL, CMA_MAX, AR_MIN, AR_MAX, TIP_CHORD_MIN,
    planform, total_weight, estimate_cd0, generate_wing, _run_vspaero, _pick,
)

N_SAMPLES = 256
OAT_NPTS = 11
ACTIVE_TOL = 0.05        # normalized margin below this counts as active
BEST_DECILE = 0.10
KEEP_PCTL = (2.0, 98.0)  # percentiles of feasible samples defining tightened bounds
FREEZE_FRAC = 0.02       # OAT drag range below this fraction of median drag -> freeze

SAMPLE_CSV = "dse_samples.csv"
OAT_CSV = "dse_oat.csv"

CONSTRAINTS = ["lift", "trim", "cma", "ar_min", "ar_max", "tip_chord"]
DERIVED = ["AR", "S_ref", "MAC", "tip_chord", "CL", "CD", "LD", "drag", "lift",
           "weight", "CM_cg", "CM_alpha_cg", "x_np", "x_cg", "feasible"]
FIELDS = VAR_NAMES + DERIVED + [f"v_{c}" for c in CONSTRAINTS]

LO = np.array([b[0] for b in BOUNDS])
HI = np.array([b[1] for b in BOUNDS])
ND = len(VAR_NAMES)

def evaluate(x):
    root_chord, taper, sweep, washout, span, dihedral, alpha = [float(v) for v in x]
    run_id = f"dse_{uuid.uuid4().hex[:8]}"
    try:
        g = planform(root_chord, taper, span)
        s_ref, ar, mac = g["s_ref"], g["ar"], g["mac"]
        weight = total_weight(span, s_ref)
        cd0 = estimate_cd0(mac, s_ref, sweep)

        vsp3_path, _ = generate_wing(run_id, span, root_chord, taper,
                                     sweep, dihedral, washout)
        r = _run_vspaero(vsp3_path, s_ref, span, mac, x_cg=0.0, cd0=cd0,
                         alpha_start=alpha - ALPHA_STEP_DEG,
                         alpha_end=alpha + ALPHA_STEP_DEG, alpha_npts=3)

        a = np.asarray(r["Alpha"], dtype=float)
        cl_m, cl_0, cl_p = [_pick(r["CL"], a, alpha + d * ALPHA_STEP_DEG) for d in (-1, 0, 1)]
        cm_m, cm_0, cm_p = [_pick(r["CMy"], a, alpha + d * ALPHA_STEP_DEG) for d in (-1, 0, 1)]
        cd_0 = _pick(r["CD"], a, alpha)

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
            "cma": (cm_alpha_cg - CMA_MAX) / abs(CMA_MAX),
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
    for i in range(done, N_SAMPLES):
        row = evaluate(pts[i])
        if row is None:
            continue
        append(SAMPLE_CSV, row)
        tag = "FEASIBLE" if row["feasible"] else "infeasible"
        print(f"  [{i + 1}/{N_SAMPLES}] {tag:10s} D={row['drag']:.4f} N  "
              f"L/D={row['LD']:.2f}  AR={row['AR']:.2f}")

def baseline(d):
    m = d["feasible"] > 0.5
    if m.sum() < 3:
        return 0.5 * (LO + HI)
    x = np.array([np.median(d[v][m]) for v in VAR_NAMES])
    return np.clip(x, LO, HI)

def oat(d):
    if count_rows(OAT_CSV) >= ND * OAT_NPTS:
        print("[oat] sweeps already on disk; skipping.")
        return
    x0 = baseline(d)
    print(f"[oat] baseline = {dict(zip(VAR_NAMES, np.round(x0, 4)))}")
    for j, name in enumerate(VAR_NAMES):
        for val in np.linspace(LO[j], HI[j], OAT_NPTS):
            x = x0.copy()
            x[j] = val
            row = evaluate(x)
            if row is None:
                continue
            append(OAT_CSV, row)
            print(f"  [{name}={val:.3f}] D={row['drag']:.4f} N  CM_cg={row['CM_cg']:+.4f}")

def oat_ranges(o, d):
    x0 = baseline(d)
    ranges = {}
    for j, name in enumerate(VAR_NAMES):
        m = np.ones(len(o["drag"]), dtype=bool)
        for k, other in enumerate(VAR_NAMES):
            if k != j:
                m &= np.isclose(o[other], x0[k], rtol=0, atol=1e-6)
        dr = o["drag"][m]
        ranges[name] = float(np.nanmax(dr) - np.nanmin(dr)) if dr.size else 0.0
    return ranges

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

def activity(d):
    m = d["feasible"] > 0.5
    drag = d["drag"][m]
    if drag.size == 0:
        return {}
    cut = np.quantile(drag, BEST_DECILE)
    best = drag <= cut
    out = {}
    for c in CONSTRAINTS:
        v = d[f"v_{c}"][m]
        act = np.abs(v) <= ACTIVE_TOL
        out[c] = {"all": float(act.mean()), "best": float(act[best].mean())}
    return out

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

def plot_surfaces(d, ranges):
    m = d["feasible"] > 0.5
    if m.sum() < 20:
        print("[plot] too few feasible samples for response surfaces")
        return
    top = sorted(VAR_NAMES, key=lambda n: -ranges.get(n, 0.0))[:3]
    pairs = [(top[0], top[1]), (top[0], top[2]), (top[1], top[2])]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for a, (xn, yn) in zip(axes, pairs):
        i, j = VAR_NAMES.index(xn), VAR_NAMES.index(yn)
        p = np.column_stack([(d[xn][m] - LO[i]) / (HI[i] - LO[i]),
                             (d[yn][m] - LO[j]) / (HI[j] - LO[j])])
        rbf = RBFInterpolator(p, d["drag"][m], kernel="thin_plate_spline", smoothing=1e-3)
        gx, gy = np.meshgrid(np.linspace(0, 1, 80), np.linspace(0, 1, 80))
        q = np.column_stack([gx.ravel(), gy.ravel()])
        z = rbf(q).reshape(gx.shape)
        z[cKDTree(p).query(q)[0].reshape(gx.shape) > 0.09] = np.nan
        cf = a.contourf(LO[i] + gx * (HI[i] - LO[i]), LO[j] + gy * (HI[j] - LO[j]),
                        z, levels=18, cmap="viridis")
        a.scatter(d[xn][m], d[yn][m], s=5, c="w", lw=0, alpha=0.5)
        a.set_xlabel(xn)
        a.set_ylabel(yn)
        fig.colorbar(cf, ax=a, label="drag [N]")
    fig.suptitle("Drag response surfaces over feasible samples", fontsize=11)
    fig.tight_layout()

def summarize(d, box, ranges, act):
    m = d["feasible"] > 0.5
    med = float(np.median(d["drag"][m])) if m.sum() else float("nan")
    x0 = baseline(d)
    freeze = {n for n in VAR_NAMES if ranges.get(n, 0.0) < FREEZE_FRAC * med}
    kept = in_box(d, box) & m

    print("\n--- FEASIBILITY ---")
    print(f"samples            : {len(d['drag'])}")
    print(f"feasible fraction  : {m.mean() * 100:.1f}%")
    print(f"after tightening   : {kept.sum() / max(in_box(d, box).sum(), 1) * 100:.1f}% "
          f"({kept.sum() / max(m.sum(), 1) * 100:.0f}% of feasible designs retained)")

    print("\n--- TIGHTENED BOUNDS ---")
    for i, n in enumerate(VAR_NAMES):
        print(f"{n:12s}: [{LO[i]:7.3f}, {HI[i]:7.3f}] -> [{box[n][0]:7.3f}, {box[n][1]:7.3f}]")

    print("\n--- OAT SENSITIVITY (drag range over full sweep) ---")
    for n, v in sorted(ranges.items(), key=lambda kv: -kv[1]):
        flag = f"  <- freeze at {x0[VAR_NAMES.index(n)]:.3f}" if n in freeze else ""
        print(f"{n:12s}: {v:.4f} N  ({v / med * 100 if med else 0:5.1f}% of median drag){flag}")

    print("\n--- CONSTRAINT ACTIVITY ---")
    for c in CONSTRAINTS:
        if c not in act:
            continue
        b = act[c]["best"]
        verdict = "binding" if b > 0.7 else "occasionally active" if b > 0.15 else "never active"
        print(f"{c:12s}: all {act[c]['all'] * 100:5.1f}%   best-decile {b * 100:5.1f}%   {verdict}")
    print()

def main():
    sample()
    d = load(SAMPLE_CSV)
    if d is None:
        print("No samples were logged.")
        return
    oat(d)
    o = load(OAT_CSV)
    ranges = oat_ranges(o, d) if o is not None else {}

    box = tighten(d)
    act = activity(d)

    summarize(d, box, ranges, act)

    plot_pairs(d, box)
    plot_surfaces(d, ranges)
    plt.show()

if __name__ == "__main__":
    main()