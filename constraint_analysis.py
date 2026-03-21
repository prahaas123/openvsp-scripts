import numpy as np
import matplotlib.pyplot as plt

# Design Requirements
V_cruise = 20.0          # Target cruise velocity (m/s)
s_to = 10.0              # Maximum allowable takeoff ground roll (m)
s_land = 30.0            # Maximum allowable landing ground roll (m)
k_L = 0.255              # Landing constant (kg/m^3)

# Parameters
g = 9.81
rho = 1.225
mu = 0.3 # Rolling friction coefficient
C_Lmax = 1.5
C_D0 = 0.02 # Parastic drag coefficient
AR = 6.0
e = 0.8 # Oswald efficiency factor
K = 1 / (np.pi * e * AR) # Induced drag factor

# Wing Loading Array
W_S = np.linspace(0, 125, 500)

# Constraint Equations
WS_landing_max = k_L * 1.0 * C_Lmax * s_land * g # Landing Constraint
TW_takeoff = (1.44 * W_S) / (rho * g * C_Lmax * s_to) + mu # Takeoff Constraint
q = 0.5 * rho * V_cruise**2
TW_cruise = (q * C_D0) / W_S + (K / q) * W_S # Cruise Constraint

# Diagram Plot
plt.figure(figsize=(10, 6))
plt.plot(W_S, TW_takeoff, label='Takeoff Constraint', color='blue', linewidth=2)
plt.plot(W_S, TW_cruise, label='Cruise Constraint', color='orange', linewidth=2)
plt.axvline(x=WS_landing_max, color='red', linestyle='--', label='Landing Constraint', linewidth=2)
valid_WS = W_S[W_S <= WS_landing_max]
valid_TW_takeoff = TW_takeoff[W_S <= WS_landing_max]
valid_TW_cruise = TW_cruise[W_S <= WS_landing_max]
bottom_boundary = np.maximum(valid_TW_takeoff, valid_TW_cruise)
plt.fill_between(valid_WS, bottom_boundary, 2.0, color='green', alpha=0.2, label='Feasible Design Space')

plt.title('Constraint Analysis Diagram (T/W vs W/S)', fontsize=14, fontweight='bold')
plt.xlabel('Wing Loading, W/S (N/m²)', fontsize=12)
plt.ylabel('Thrust-to-Weight Ratio, T/W', fontsize=12)
plt.ylim(0, 2.0)
plt.xlim(0, max(W_S))
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper left', framealpha=1)

plt.tight_layout()
plt.show()