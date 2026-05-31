import cv2
import numpy as np
from scipy.interpolate import splprep, splev

# INPUTS
INPUT_IMAGE  = "Tip.png"
OUTPUT_NAME  = "tip"
NUM_POINTS   = 200             # total output points (split across upper/lower)
THRESHOLD    = 80              # 0–255, lower = pick up lighter lines

def load_and_threshold(path, thresh):
    img  = cv2.imread(path)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    _, bw = cv2.threshold(cv2.bitwise_not(gray), thresh, 255, cv2.THRESH_BINARY)
    return bw

def largest_contour(bw):
    contours, _ = cv2.findContours(bw, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    return max(contours, key=lambda c: cv2.arcLength(c, closed=True))

def contour_to_xy(contour):
    pts = contour[:, 0, :]
    return pts[:, 0].astype(float), -pts[:, 1].astype(float)  # flip y

def resample_spline(x, y, n):
    tck, _ = splprep([x, y], s=0, per=True, k=3)
    xr, yr = splev(np.linspace(0, 1, n, endpoint=False), tck)
    return np.array(xr), np.array(yr)

def normalize(x, y):
    i_le = np.argmin(x)
    x = x - x[i_le];  y = y - y[i_le]
    i_te  = np.argmax(x)
    angle = np.arctan2(y[i_te], x[i_te])
    c, s  = np.cos(-angle), np.sin(-angle)
    x, y  = c*x - s*y, s*x + c*y
    return x / x.max(), y / x.max()

def split_and_order(x, y, n):
    i_le = np.argmin(x)
    i_te = np.argmax(x)

    # Two paths between LE and TE
    if i_le < i_te:
        seg_a = np.column_stack((x[i_le:i_te+1], y[i_le:i_te+1]))
        seg_b = np.column_stack((np.r_[x[i_te:], x[:i_le+1]],
                                 np.r_[y[i_te:], y[:i_le+1]]))
    else:
        seg_a = np.column_stack((x[i_te:i_le+1], y[i_te:i_le+1]))
        seg_b = np.column_stack((np.r_[x[i_le:], x[:i_te+1]],
                                 np.r_[y[i_le:], y[:i_te+1]]))

    # Upper surface has positive mean y
    upper, lower = (seg_a, seg_b) if seg_a[:,1].mean() > seg_b[:,1].mean() else (seg_b, seg_a)

    # Ensure upper runs LE→TE (x increasing), lower runs LE→TE too
    if upper[0,0] > upper[-1,0]: upper = upper[::-1]
    if lower[0,0] > lower[-1,0]: lower = lower[::-1]

    def reseg(seg, k):
        t0 = np.linspace(0, 1, len(seg))
        t  = np.linspace(0, 1, k)
        return np.column_stack((np.interp(t, t0, seg[:,0]),
                                np.interp(t, t0, seg[:,1])))

    half = n // 2
    upper_rs = reseg(upper[::-1], half)   # TE → LE
    lower_rs  = reseg(lower,       half)  # LE → TE
    return np.vstack([upper_rs, lower_rs[1:]])   # drop duplicate LE

def write_dat(coords, name):
    filename = name + ".dat"
    with open(filename, "w") as f:
        f.write(name + "\n")
        for xi, yi in coords:
            f.write(f"{xi:.6f}  {yi:.6f}\n")
    print(f"Saved {filename}  ({len(coords)} points)")

# MAIN
bw      = load_and_threshold(INPUT_IMAGE, THRESHOLD)
contour = largest_contour(bw)
x, y    = contour_to_xy(contour)
x, y    = resample_spline(x, y, 2000)
x, y    = normalize(x, y)
coords  = split_and_order(x, y, NUM_POINTS)
write_dat(coords, OUTPUT_NAME)