import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

# Aircraft constants
mass = 1.35                  # kg
I_y = 0.35                  # kg*m^2 (estimated realistic value)
Cmq = -0.5                  # Pitch damping derivative (estimated)

rho = 1.225                 # kg/m^3
g = 9.81                    # m/s^2

aero_file = "aero_full.csv"
df = pd.read_csv(aero_file)
points = np.column_stack((df["Velocity"], df["Alpha_deg"]))
f_CL = LinearNDInterpolator(points, df["CL"])
f_CD = LinearNDInterpolator(points, df["CD"])
f_Cm = LinearNDInterpolator(points, df["Cm"])
V_min_db = df["Velocity"].min()
V_max_db = df["Velocity"].max()
alpha_min_db = df["Alpha_deg"].min()
alpha_max_db = df["Alpha_deg"].max()

stability_file = "stability.csv"
df_stab = pd.read_csv(stability_file)

# Reference geometry the VSPAERO coefficients were normalised with
S_wing = df_stab["Sref"].iloc[0]    # m^2
c_bar = df_stab["cref"].iloc[0]     # m

def make_interp(x, y):
    # interp1d needs >= 2 points; a single-velocity stability.csv is treated as constant
    if len(x) == 1:
        val = float(y.iloc[0])
        return lambda V: val
    return interp1d(x, y, kind='linear', fill_value='extrapolate')

f_Cmq = make_interp(df_stab["Velocity"], df_stab["Cm_q"])

# Finding trim conditions
print("Searching for approximate trim condition...")
best_error = 1e9
trim_index = 0

for i in range(len(df)):
    V = df["Velocity"].iloc[i]
    alpha_deg = df["Alpha_deg"].iloc[i]
    CL = df["CL"].iloc[i]
    Cm = df["Cm"].iloc[i]
    q_dyn = 0.5 * rho * V**2
    Lift = q_dyn * S_wing * CL
    lift_error = abs(Lift - mass * g)
    moment_error = abs(Cm)
    total_error = lift_error + 10 * moment_error

    if total_error < best_error:
        best_error = total_error
        trim_index = i

# Extract trim state
V_trim = df["Velocity"].iloc[trim_index]
alpha_trim_deg = df["Alpha_deg"].iloc[trim_index]
alpha_trim = np.deg2rad(alpha_trim_deg)
CL_trim = df["CL"].iloc[trim_index]
CD_trim = df["CD"].iloc[trim_index]
q_dyn_trim = 0.5 * rho * V_trim**2
Lift_trim = q_dyn_trim * S_wing * CL_trim
Drag_trim = q_dyn_trim * S_wing * CD_trim
T_trim = Drag_trim # Assume trim thrust equals trim drag

print(f"Trim velocity: {V_trim:.2f} m/s")
print(f"Trim alpha: {alpha_trim_deg:.2f} deg")

# Initial disturbance from trim condition
alpha0 = alpha_trim + np.deg2rad(3.0)   # 3 deg pitch disturbance
q0 = np.deg2rad(5.0)                    # initial pitch rate
theta0 = alpha0                         # assume initial gamma = 0
V0 = V_trim
x0 = [V0, alpha0, q0, theta0]

def longitudinal_dynamics(t, x):
    V, alpha, q, theta = x
    V = max(V, 1.0)
    alpha_deg = np.rad2deg(alpha)
    # Coefficients are nearly Reynolds-independent, so clamp V to the swept range for lookup only
    V_lookup = np.clip(V, V_min_db, V_max_db)
    CL = f_CL(V_lookup, alpha_deg)
    CD = f_CD(V_lookup, alpha_deg)
    Cm = f_Cm(V_lookup, alpha_deg)

    # Check for NaN
    if np.isnan(CL) or np.isnan(CD) or np.isnan(Cm):

        print("\nWARNING:")
        print(f"State outside aerodynamic database at t = {t:.2f}s")
        print(f"V = {V:.2f} m/s")
        print(f"alpha = {alpha_deg:.2f} deg")

        return [0, 0, 0, 0]

    Q = 0.5 * rho * V**2
    Cmq = f_Cmq(V)
    Cm_total = Cm + Cmq * (q * c_bar / (2 * V))
    Lift = Q * S_wing * CL
    Drag = Q * S_wing * CD
    M_aero = Q * S_wing * c_bar * Cm_total

    gamma = theta - alpha

    # Equations of motion
    V_dot = (T_trim - Drag) / mass - g * np.sin(gamma)
    alpha_dot = q - ((Lift - mass * g * np.cos(gamma)) / (mass * V))
    q_dot = M_aero / I_y
    theta_dot = q

    return [V_dot, alpha_dot, q_dot, theta_dot]

# Stop the integration once alpha leaves the aerodynamic database
def out_of_database(t, x):
    alpha_deg = np.rad2deg(x[1])
    return min(alpha_deg - alpha_min_db, alpha_max_db - alpha_deg)
out_of_database.terminal = True

# Simulation

t_start = 0
t_end = 30

t_eval = np.linspace(t_start, t_end, 3000)

print("Running simulation...")

solution = solve_ivp(
    longitudinal_dynamics,
    [t_start, t_end],
    x0,
    method='RK45',
    t_eval=t_eval,
    events=out_of_database,
    rtol=1e-6,
    atol=1e-8
)

if solution.status == 1:
    print(f"\nWARNING: alpha left the aerodynamic database "
          f"({alpha_min_db}..{alpha_max_db} deg) at t = {solution.t_events[0][0]:.2f}s; simulation stopped.")

t = solution.t
V = solution.y[0]
alpha = np.rad2deg(solution.y[1])
q = np.rad2deg(solution.y[2])
theta = np.rad2deg(solution.y[3])
gamma = theta - alpha

# Plot results

plt.figure(figsize=(12, 10))

plt.subplot(4, 1, 1)
plt.plot(t, V)
plt.ylabel("Velocity (m/s)")
plt.grid(True)
plt.title("Longitudinal Stability Response")

plt.subplot(4, 1, 2)
plt.plot(t, alpha)
plt.ylabel("Alpha (deg)")
plt.grid(True)

plt.subplot(4, 1, 3)
plt.plot(t, q)
plt.ylabel("Pitch Rate (deg/s)")
plt.grid(True)

plt.subplot(4, 1, 4)
plt.plot(t, theta, label='Theta')
plt.plot(t, gamma, label='Gamma')
plt.xlabel("Time (s)")
plt.ylabel("Angle (deg)")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

# Print final stability summary

final_alpha = alpha[-1]
final_q = q[-1]

print("\n==============================")
print("STABILITY SUMMARY")
print("==============================")

print(f"Final alpha: {final_alpha:.3f} deg")
print(f"Final pitch rate: {final_q:.3f} deg/s")

if abs(final_q) < 0.5 and abs(final_alpha - alpha_trim_deg) < 1.0:
    print("Aircraft appears dynamically stable.")
else:
    print("Aircraft may be unstable or weakly damped.")