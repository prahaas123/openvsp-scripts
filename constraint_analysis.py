import numpy as np
import matplotlib.pyplot as plt

# Design Requirements
V_cruise = 10.0          # Target cruise velocity (m/s)
s_to = 5.0               # Maximum allowable takeoff ground roll (m)
V_stall_max = 8.0        # Maximum allowable stall speed (m/s)

# Parameters
g = 9.81
rho = 1.225
mu = 0.05 # Rolling friction coefficient (~0.04 pavement, ~0.08 grass)
C_Lmax = 1.5
C_L_ground = 0.3 # Lift coefficient at ground-roll attitude
C_D0 = 0.04 # Parasitic drag coefficient, whole aircraft (wing alone is only ~0.013)
AR = 5.0
e = 0.8 # Oswald efficiency factor
K = 1 / (np.pi * e * AR) # Induced drag factor
k_LO = 1.2               # Liftoff speed / stall speed

# Wing Loading Array
W_S = np.linspace(1, 100, 500)
TW_plot_max = 2.0

# Constraint Equations
CD_ground = C_D0 + K * C_L_ground**2
q_to = 0.49 * k_LO**2 / C_Lmax # q / (W/S) during the takeoff roll
TW_takeoff = (k_LO**2 * W_S) / (rho * g * C_Lmax * s_to) + mu + (CD_ground - mu * C_L_ground) * q_to # Takeoff Constraint
WS_stall_max = 0.5 * rho * V_stall_max**2 * C_Lmax # Stall Speed Constraint
q = 0.5 * rho * V_cruise**2
TW_cruise = (q * C_D0) / W_S + (K / q) * W_S # Cruise Constraint

# Diagram Plot
plt.figure(figsize=(10, 6))
plt.plot(W_S, TW_takeoff, label='Takeoff Constraint', color='blue', linewidth=2)
plt.plot(W_S, TW_cruise, label='Cruise Constraint', color='orange', linewidth=2)
plt.axvline(x=WS_stall_max, color='black', linestyle='--', label='Stall Speed Constraint', linewidth=2)

valid = W_S <= WS_stall_max
bottom_boundary = np.maximum(TW_takeoff[valid], TW_cruise[valid])
plt.fill_between(W_S[valid], bottom_boundary, TW_plot_max, where=bottom_boundary <= TW_plot_max,
                 color='green', alpha=0.2, label='Feasible Design Space')

plt.title('Constraint Analysis Diagram (T/W vs W/S)', fontsize=14, fontweight='bold')
plt.xlabel('Wing Loading, W/S (N/m²)', fontsize=12)
plt.ylabel('Thrust-to-Weight Ratio, T/W', fontsize=12)
plt.ylim(0, TW_plot_max)
plt.xlim(0, max(W_S))
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper left', framealpha=1)

plt.tight_layout()
plt.show()
