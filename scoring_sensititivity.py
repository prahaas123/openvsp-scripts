import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

# ==========================================
# 1. BASELINE PARAMETERS
# ==========================================
W_empty_base = 2.0     # Empty weight (lbs)
W_payload_base = 1.5   # Payload weight (lbs)
S_ft_base = 4.0        # Wingspan (ft)
V_LOF = 12.5           # Lift-off velocity (m/s)
rho = 1.225            # Air density (kg/m^3)
mu = 0.05              # Rolling friction
C_L = 0.8              # Lift coefficient
C_D = 0.1              # Drag coefficient
c_root = 0.9           # Wing root chord (m)

def evaluate_score(k_We=1.0, k_Wp=1.0, k_L=1.0, k_D=1.0, k_S=1.0):
    """
    Calculates the final flight score by first solving the takeoff physics
    integration, and then feeding that distance into the scoring equation.
    """
    # Adjust parameters based on multipliers
    W_e_lbs = W_empty_base * k_We
    W_p_lbs = W_payload_base * k_Wp
    S_ft = S_ft_base * k_S
    
    # Convert to metric for Physics Integration
    W_total_N = (W_e_lbs + W_p_lbs) * 4.44822  # 1 lb = 4.448 N
    mass_kg = W_total_N / 9.81
    
    b_m = S_ft * 0.3048                        # Wingspan in meters
    Wing_Area = 0.5 * b_m * c_root             # Delta wing area (m^2)
    
    # ==========================================
    # 2. PHYSICS INTEGRATION (Takeoff Distance)
    # ==========================================
    def integrand(V):
        T = 30.0 - 0.5 * V
        L = k_L * 0.5 * rho * V**2 * Wing_Area * C_L
        D = k_D * 0.5 * rho * V**2 * Wing_Area * C_D
        
        # Normal force cannot be negative (plane is off the ground)
        N_force = max(0.0, W_total_N - L)
        net_force = T - D - mu * N_force
        
        # Prevent division by zero if acceleration stalls
        if net_force <= 0.01:
            return 1e6
        return V / net_force
        
    # Integrate to find distance in meters, then convert to feet
    dist_m, _ = quad(integrand, 0, V_LOF)
    dist_m = mass_kg * dist_m
    dist_ft = dist_m / 0.3048
    
    # ==========================================
    # 3. SCORING EQUATION
    # ==========================================
    M = 11.0 / ((W_e_lbs - 1.0)**4 + 8.9)
    
    # Linear approximation of Takeoff Bonus: capped at 0 minimum
    B_takeoff = max(0.0, 20.0 - 0.2 * dist_ft)
    Z = B_takeoff - S_ft**1.5
    
    # Final Flight Score calculation
    FS = 3.0 * W_p_lbs * M + Z
    return FS

# ==========================================
# 4. RUN SENSITIVITY SWEEP
# ==========================================
multipliers = np.linspace(0.8, 1.2, 50)
pct_changes = (multipliers - 1) * 100

base_score = evaluate_score()
print(f"Baseline Score: {base_score:.2f} points")

scores = {'We': [], 'Wp': [], 'L': [], 'D': [], 'S': []}

for k in multipliers:
    scores['We'].append(evaluate_score(k_We=k))
    scores['Wp'].append(evaluate_score(k_Wp=k))
    scores['L'].append(evaluate_score(k_L=k))
    scores['D'].append(evaluate_score(k_D=k))
    scores['S'].append(evaluate_score(k_S=k))

# Convert to percentage changes
plt.figure(figsize=(10,6))
plt.plot(pct_changes, (np.array(scores['We']) - base_score)/base_score * 100, label='Empty Weight ($W_{Empty}$)', color='red', linewidth=2)
plt.plot(pct_changes, (np.array(scores['Wp']) - base_score)/base_score * 100, label='Payload Weight ($W_{Payload}$)', color='green', linewidth=2)
plt.plot(pct_changes, (np.array(scores['L']) - base_score)/base_score * 100, label='Lift', color='blue', linewidth=2)
plt.plot(pct_changes, (np.array(scores['D']) - base_score)/base_score * 100, label='Drag', color='orange', linewidth=2)
plt.plot(pct_changes, (np.array(scores['S']) - base_score)/base_score * 100, label='Wingspan ($S$)', color='purple', linewidth=2)

plt.axhline(0, color='k', linestyle='--', alpha=0.5)
plt.axvline(0, color='k', linestyle='--', alpha=0.5)
plt.xlabel('Change in Parameter (%)', fontsize=12)
plt.ylabel('Change in Flight Score (%)', fontsize=12)
plt.title('Integrated Sensitivity Analysis: Flight Score vs. Design Parameters', fontsize=14)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()