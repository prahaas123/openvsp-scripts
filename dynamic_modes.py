
import numpy as np
import matplotlib.pyplot as plt
import cmath as cm
import pandas as pd

plt.rcParams.update({'font.size': 14})

# Aircraft mass & inertia properties
mass   = 1.35            # kg
Ix     = 0.15            # kg*m^2
Iy     = 0.35            # kg*m^2
Iz     = 0.45            # kg*m^2
Ixz    = 0.00            # kg*m^2

# Assumed stability derivatives (per rad)
CLq = 6.000     # CL_q
CDq = 0.000     # CD_q

# Constants
g0      = 9.81    # m/s^2
density = 1.225   # kg/m^3

# Reference Geometry (the values the VSPAERO derivatives were normalised with)
df_stab = pd.read_csv('stability.csv')
df_aero = pd.read_csv('aero_full.csv')
S_ref = df_stab['Sref'].iloc[0]   # m^2 (Reference Area)
c_ref = df_stab['cref'].iloc[0]   # m   (Reference Chord)
b_ref = df_stab['bref'].iloc[0]   # m   (Reference Span)

# Aero polar grid: rows = velocity, columns = alpha
CL_grid = df_aero.pivot(index='Velocity', columns='Alpha_deg', values='CL')
CD_grid = df_aero.pivot(index='Velocity', columns='Alpha_deg', values='CD')
CM_grid = df_aero.pivot(index='Velocity', columns='Alpha_deg', values='Cm')
vel_grid   = CL_grid.index.values.astype(float)
alpha_grid = CL_grid.columns.values.astype(float)   # deg

def polar_at(grid, V):
    # Coefficient vs alpha at velocity V (linear in velocity between sweep points)
    return np.array([np.interp(V, vel_grid, grid[a].values) for a in alpha_grid])

def interp_extrap(x, xp, fp):
    # Linear interpolation, extrapolating with the end slopes outside xp
    if x < xp[0]:
        return fp[0] + (fp[1] - fp[0]) / (xp[1] - xp[0]) * (x - xp[0])
    if x > xp[-1]:
        return fp[-1] + (fp[-1] - fp[-2]) / (xp[-1] - xp[-2]) * (x - xp[-1])
    return np.interp(x, xp, fp)

# Trim Conditions
Vtrim = float(df_stab['Velocity'].iloc[0])   # m/s, initial guess
for _ in range(50):
    CL_polar = polar_at(CL_grid, Vtrim)
    CM_polar = polar_at(CM_grid, Vtrim)
    # Cm decreases with alpha, reverse so the interpolation x-values increase
    alpha_trim_deg = float(interp_extrap(0.0, CM_polar[::-1], alpha_grid[::-1]))
    CL0 = float(interp_extrap(alpha_trim_deg, alpha_grid, CL_polar))
    V_new = np.sqrt(2.0 * mass * g0 / (density * S_ref * CL0))
    converged = abs(V_new - Vtrim) < 1e-6
    Vtrim = V_new
    if converged:
        break
Q = 0.5 * density * Vtrim**2
theta0 = alpha_trim_deg * np.pi / 180.0   # gamma = 0, so theta = alpha

if not (vel_grid[0] <= Vtrim <= vel_grid[-1]):
    print(f'WARNING: trim speed {Vtrim:.2f} m/s is outside the aero_full.csv velocity range '
          f'({vel_grid[0]:.0f} to {vel_grid[-1]:.0f} m/s)')

# Stability derivatives in stability.csv are assumed to apply at Vtrim
vel_data = np.full(len(df_stab), Vtrim)

# Base Coefficients at Trim
CL_polar = polar_at(CL_grid, Vtrim)
CD_polar = polar_at(CD_grid, Vtrim)
CM_polar = polar_at(CM_grid, Vtrim)
CD0 = float(interp_extrap(alpha_trim_deg, alpha_grid, CD_polar))
CM0 = float(interp_extrap(alpha_trim_deg, alpha_grid, CM_polar))

# AERODYNAMIC COEFFICIENTS
def get_deriv(col_name):
    return np.interp(Vtrim, vel_data, df_stab[col_name].values)

def speed_deriv(grid):
    # Nondimensional speed derivative C_u = V * dC/dV at fixed trim alpha
    C_vs_V = np.array([interp_extrap(alpha_trim_deg, alpha_grid, grid.loc[v].values)
                       for v in vel_grid])
    return Vtrim * np.interp(Vtrim, vel_grid, np.gradient(C_vs_V, vel_grid))

# Longitudinal Derivatives (per rad)
CLa = get_deriv('CL_alpha')
CMa = get_deriv('Cm_alpha')
CDa = float(interp_extrap(alpha_trim_deg, alpha_grid,
                          np.gradient(CD_polar, np.deg2rad(alpha_grid))))
CLu = speed_deriv(CL_grid)
CDu = speed_deriv(CD_grid)
CMu = speed_deriv(CM_grid)
CMq = get_deriv('Cm_q')

print(f'Trim: V = {Vtrim:.2f} m/s, alpha = {alpha_trim_deg:.2f} deg')
print(f'CL0 = {CL0:.4f}, CD0 = {CD0:.4f}, CM0 = {CM0:.4f}')
print(f'CLa = {CLa:.4f}, CDa = {CDa:.4f}, CMa = {CMa:.4f} (per rad)')
print(f'CLu = {CLu:.4f}, CDu = {CDu:.4f}, CMu = {CMu:.4f}')

# Lateral Derivatives
CYb = get_deriv('CY_beta')
Clb = get_deriv('Cl_beta')
Cnb = get_deriv('Cn_beta')
CYp = get_deriv('CY_p')
Clp = get_deriv('Cl_p')
Cnp = get_deriv('Cn_p')
CYr = get_deriv('CY_r')
Clr = get_deriv('Cl_r')
Cnr = get_deriv('Cn_r')

# DIMENSIONAL DERIVATIVES
# Longitudinal Dimensional Derivatives
# X-force
XU = -(Q * S_ref * (2.0 * CD0 + CDu)) / (mass * Vtrim)
XW =  (Q * S_ref * (CL0 - CDa)) / (mass * Vtrim)
XQ = -(Q * S_ref * c_ref * CDq) / (2.0 * mass * Vtrim)

# Z-force
ZU = -(Q * S_ref * (2.0 * CL0 + CLu)) / (mass * Vtrim)
ZW = -(Q * S_ref * (CLa + CD0)) / (mass * Vtrim)
ZQ = -(Q * S_ref * c_ref * CLq) / (2.0 * mass * Vtrim)
ZWDOT = 0.0

# Pitching Moment
MU =  (Q * S_ref * c_ref * CMu) / (Iy * Vtrim)
MW =  (Q * S_ref * c_ref * CMa) / (Iy * Vtrim)
MQ =  (Q * S_ref * c_ref**2 * CMq) / (2.0 * Iy * Vtrim)
MWDOT = 0.0

# Lateral Dimensional Derivatives
# Y-force
YV =  (Q * S_ref * CYb) / (mass * Vtrim)
YP =  (Q * S_ref * b_ref * CYp) / (2.0 * mass * Vtrim)
YR =  (Q * S_ref * b_ref * CYr) / (2.0 * mass * Vtrim)

# Rolling Moment
LV =  (Q * S_ref * b_ref * Clb) / (Ix * Vtrim)
LP =  (Q * S_ref * b_ref**2 * Clp) / (2.0 * Ix * Vtrim)
LR =  (Q * S_ref * b_ref**2 * Clr) / (2.0 * Ix * Vtrim)

# Yawing Moment
NV =  (Q * S_ref * b_ref * Cnb) / (Iz * Vtrim)
NP =  (Q * S_ref * b_ref**2 * Cnp) / (2.0 * Iz * Vtrim)
NR =  (Q * S_ref * b_ref**2 * Cnr) / (2.0 * Iz * Vtrim)


# STATE SPACE MATRICES
# Longitudinal A-Matrix (States: u, w, q, theta)
first_row  = [XU, XW, 0, -g0 * np.cos(theta0)]
second_row = [ZU, ZW, Vtrim, -g0 * np.sin(theta0)]
third_row  = [MU + MWDOT * ZU, MW + MWDOT * ZW, MQ + Vtrim * MWDOT, -MWDOT * g0 * np.sin(theta0)]
fourth_row = [0., 0., 1., 0.]
ALON = [first_row, second_row, third_row, fourth_row]

# Approximation to Phugoid
Aph = [[XU, -g0], [-ZU / Vtrim, 0]]
# Approximation to Short Period
Asp = [[ZW / (1 - ZWDOT), (Vtrim + ZQ) / (1 - ZWDOT)],
       [MW + (MWDOT * ZW) / (1 - ZWDOT), MQ + MWDOT * (Vtrim + ZQ) / (1 - ZWDOT)]]

# Create Lateral A-Matrix
first_row  = [YV, YP, g0 * np.cos(theta0), YR - Vtrim]
second_row = [LV, LP, 0., LR]
third_row  = [0., 1., 0., 0.]
fourth_row = [NV, NP, 0., NR]
ALAT = [first_row, second_row, third_row, fourth_row]

# Approximation to Dutch Roll
a = 1.0
b = -((LP * NR + Vtrim * NV - LR * NP) / (LP + NR)
      + Vtrim * (LV * NP - LP * NV) / (LP + NR)**2)
c = Vtrim * (LP * NV - LV * NP) / (LP + NR)
dutch_roll1 = -b / (2 * a) + 0.5 * cm.sqrt(b**2 - 4 * a * c)
dutch_roll2 = -b / (2 * a) - 0.5 * cm.sqrt(b**2 - 4 * a * c)

# Approximation to Roll Mode
ix = Ixz / Ix
iz = Ixz / Iz
roll_mode = (LP + ix * NP) / (1 - ix * iz)

# Approximation to Spiral
spiral_mode = NR - LR * NV / LV

# COMPUTE EIGENVALUES
LONeigenvalues, LONeigenvectors = np.linalg.eig(ALON)
LATeigenvalues, LATeigenvectors = np.linalg.eig(ALAT)
PHeigenvalues,  PHeigenvectors  = np.linalg.eig(Aph)
SPeigenvalues,  PHeigenvectors  = np.linalg.eig(Asp)

# PLOT AND PRINT EVERYTHING
plt.figure(figsize=(10, 6))
plt.axhline(0, color='black', linewidth=1.2, alpha=0.8)
plt.axvline(0, color='black', linewidth=1.2, alpha=0.8)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xlabel('Real Axis')
plt.ylabel('Imaginary Axis')
plt.title('Aircraft Dynamic Modes Root Locus')

def fmt(z):
    # Eigenvalue as a string to two decimals
    z = complex(z)
    if z.imag == 0:
        return f'{z.real:.2f}'
    sign = '+' if z.imag > 0 else '-'
    return f'{z.real:.2f} {sign} {abs(z.imag):.2f}j'

print('Longitudinal Modes:')
for i in range(len(LONeigenvalues)):
    print(fmt(LONeigenvalues[i]))
    # Only label the first point to avoid duplicate legend entries
    lbl = 'Longitudinal Modes' if i == 0 else None
    plt.plot(np.real(LONeigenvalues[i]), np.imag(LONeigenvalues[i]), 'bx', markersize=8, markeredgewidth=2, label=lbl)

print('Short Period Approx:')
for i in range(len(SPeigenvalues)):
    print(fmt(SPeigenvalues[i]))
    lbl = 'Short Period Approx' if i == 0 else None
    plt.plot(np.real(SPeigenvalues[i]), np.imag(SPeigenvalues[i]), 'rx', markersize=8, markeredgewidth=2, label=lbl)

print('Phugoid Approx:')
for i in range(len(PHeigenvalues)):
    print(fmt(PHeigenvalues[i]))
    lbl = 'Phugoid Approx' if i == 0 else None
    plt.plot(np.real(PHeigenvalues[i]), np.imag(PHeigenvalues[i]), 'gx', markersize=8, markeredgewidth=2, label=lbl)

print('Lateral Modes:')
for i in range(len(LATeigenvalues)):
    print(fmt(LATeigenvalues[i]))
    lbl = 'Lateral Modes' if i == 0 else None
    plt.plot(np.real(LATeigenvalues[i]), np.imag(LATeigenvalues[i]), 'bo', markersize=8, label=lbl)

print('Dutch Roll Approx:')
print(fmt(dutch_roll1))
plt.plot(np.real(dutch_roll1), np.imag(dutch_roll1), 'ro', markersize=8, label='Dutch Roll Approx')
print(fmt(dutch_roll2))
plt.plot(np.real(dutch_roll2), np.imag(dutch_roll2), 'ro', markersize=8)

print('Roll Mode Approx:')
print(fmt(roll_mode))
plt.plot(np.real(roll_mode), np.imag(roll_mode), 'go', markersize=8, label='Roll Approx')

print('Spiral Mode Approx:')
print(fmt(spiral_mode))
plt.plot(np.real(spiral_mode), np.imag(spiral_mode), 'mo', markersize=8, label='Spiral Approx')

plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left", borderaxespad=0.)
plt.tight_layout()
plt.show()