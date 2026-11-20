import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

# Aircraft inputs
mass = 1.35                  # kg
I_x = 0.15                  # Roll inertia (kg*m^2)
I_z = 0.45                  # Yaw inertia (kg*m^2)
I_xz = 0.0                  # Cross inertia (usually small for symmetric UAV)

# Initial disturbance from trim condition
beta0_deg = 5.0
beta0 = np.deg2rad(beta0_deg)
p0 = np.deg2rad(10.0)       # roll rate disturbance
r0 = 0.0                    # yaw rate disturbance
phi0 = 0.0                  # roll angle
psi0 = 0.0                  # heading angle

rho = 1.225                 # kg/m^3
g = 9.81                    # m/s^2

stability_file = "stability.csv"
df = pd.read_csv(stability_file)

# Reference geometry the VSPAERO derivatives were normalised with
S_wing = df["Sref"].iloc[0] # m^2
b_span = df["bref"].iloc[0] # m
def make_interp(x, y):
    # interp1d needs >= 2 points; a single-velocity stability.csv is treated as constant
    if len(x) == 1:
        val = float(y.iloc[0])
        return lambda V: val
    return interp1d(x, y, fill_value='extrapolate')

f_CY_beta = make_interp(df["Velocity"], df["CY_beta"])
f_CY_p    = make_interp(df["Velocity"], df["CY_p"])
f_CY_r    = make_interp(df["Velocity"], df["CY_r"])
f_Cl_beta = make_interp(df["Velocity"], df["Cl_beta"])
f_Cl_p    = make_interp(df["Velocity"], df["Cl_p"])
f_Cl_r    = make_interp(df["Velocity"], df["Cl_r"])
f_Cn_beta = make_interp(df["Velocity"], df["Cn_beta"])
f_Cn_p    = make_interp(df["Velocity"], df["Cn_p"])
f_Cn_r    = make_interp(df["Velocity"], df["Cn_r"])

# Trim conditions
V_trim = 20.0
Q_trim = 0.5 * rho * V_trim**2
v0 = V_trim * np.sin(beta0)
x0 = [v0, p0, r0, phi0, psi0]
print(f"Simulation velocity: {V_trim:.2f} m/s")

Inertia = np.array([
    [I_x, -I_xz],
    [-I_xz, I_z]
])

def lateral_dynamics(t, x):
    v, p, r, phi, psi = x
    beta = np.arcsin(np.clip(v / V_trim, -0.99, 0.99))
    p_hat = p * b_span / (2 * V_trim)
    r_hat = r * b_span / (2 * V_trim)
    
    CY_beta = f_CY_beta(V_trim)
    CY_p    = f_CY_p(V_trim)
    CY_r    = f_CY_r(V_trim)

    Cl_beta = f_Cl_beta(V_trim)
    Cl_p    = f_Cl_p(V_trim)
    Cl_r    = f_Cl_r(V_trim)

    Cn_beta = f_Cn_beta(V_trim)
    Cn_p    = f_Cn_p(V_trim)
    Cn_r    = f_Cn_r(V_trim)
    
    CY = (CY_beta * beta + CY_p * p_hat + CY_r * r_hat)
    Cl = ( Cl_beta * beta + Cl_p * p_hat + Cl_r * r_hat)
    Cn = (Cn_beta * beta + Cn_p * p_hat + Cn_r * r_hat)
    Y_force = Q_trim * S_wing * CY
    L_moment = Q_trim * S_wing * b_span * Cl
    N_moment = Q_trim * S_wing * b_span * Cn

    v_dot = (Y_force / mass - r * V_trim + g * phi)
    moments = np.array([L_moment, N_moment])
    ang_accel = np.linalg.solve(Inertia, moments)
    p_dot = ang_accel[0]
    r_dot = ang_accel[1]
    phi_dot = p
    psi_dot = r
    return [v_dot, p_dot, r_dot, phi_dot, psi_dot]

# Simulation
t_start = 0
t_end = 20
t_eval = np.linspace(t_start, t_end, 3000)
print("Running lateral-directional simulation...")
solution = solve_ivp(
    lateral_dynamics,
    [t_start, t_end],
    x0,
    method='RK45',
    t_eval=t_eval,
    rtol=1e-7,
    atol=1e-9
)

# Simulation results
t = solution.t
v = solution.y[0]
p = np.rad2deg(solution.y[1])
r = np.rad2deg(solution.y[2])
phi = np.rad2deg(solution.y[3])
psi = np.rad2deg(solution.y[4])
beta = np.rad2deg(np.arcsin(np.clip(v / V_trim, -0.99, 0.99)))

# Plot results
plt.figure(figsize=(12, 12))
plt.subplot(5, 1, 1)
plt.plot(t, beta)
plt.ylabel("Beta (deg)")
plt.grid(True)
plt.title("Lateral-Directional Stability Response")

plt.subplot(5, 1, 2)
plt.plot(t, p)
plt.ylabel("Roll Rate p (deg/s)")
plt.grid(True)

plt.subplot(5, 1, 3)
plt.plot(t, r)
plt.ylabel("Yaw Rate r (deg/s)")
plt.grid(True)

plt.subplot(5, 1, 4)
plt.plot(t, phi)
plt.ylabel("Roll Angle φ (deg)")
plt.grid(True)

plt.subplot(5, 1, 5)
plt.plot(t, psi)
plt.ylabel("Heading ψ (deg)")
plt.xlabel("Time (s)")
plt.grid(True)
plt.tight_layout()
plt.show()

# Print final stability assessment
final_beta = beta[-1]
final_p = p[-1]
final_r = r[-1]

print("\n===================================")
print("LATERAL-DIRECTIONAL STABILITY")
print("===================================")
print(f"Final beta: {final_beta:.3f} deg")
print(f"Final roll rate: {final_p:.3f} deg/s")
print(f"Final yaw rate: {final_r:.3f} deg/s")

# Dutch roll / spiral qualitative check
stable = (
    abs(final_beta) < 0.5 and
    abs(final_p) < 0.5 and
    abs(final_r) < 0.5
)

if stable:
    print("Aircraft appears laterally stable.")
else:
    print("Aircraft may be weakly damped or unstable.")

print("\n===================================")
print("DERIVATIVE CHECK")
print("===================================")
print(f"Cl_beta = {float(f_Cl_beta(V_trim)):.3f}")
print(f"Cn_beta = {float(f_Cn_beta(V_trim)):.3f}")
print(f"Cl_p    = {float(f_Cl_p(V_trim)):.3f}")
print(f"Cn_r    = {float(f_Cn_r(V_trim)):.3f}")
print("\nExpected stable signs:")
print("Cl_beta < 0  -> Dihedral roll stability")
print("Cn_beta > 0  -> Directional weathercock stability")
print("Cl_p    < 0  -> Roll damping")
print("Cn_r    < 0  -> Yaw damping")