import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_discrete_are # Import DARE solver
import matplotlib.pyplot as plt

# Set default font size for plots
plt.rcParams['font.size'] = 14 

# 1. Nonlinear System Dynamics
def nonlinear_system(x, u):
    # Scalar system dynamics: dx/dt = -x + x^3 + u
    return -x + x**3 + u

# 2. Linear Controller Parameters (derived above)
Ad = 0.905
Bd = 0.095
K = 4.263
Ts = 0.1
x_eq = 0.0 # Equilibrium state

# 3. Simulation Setup
t_start = 0
t_end = 1.0 # Shorter simulation time might be enough
n_steps = int(t_end / Ts)

# Initial Conditions to compare
x0_small = 2
x0_large = 2.4

# --- DLQR Controller Design ---
Qd = 10.0  # State weight
Rd = 0.1  # Control weight

# Solve Discrete Algebraic Riccati Equation (DARE)
# Need to pass arguments as 2D arrays for the solver
A_dare = np.array([[Ad]])
B_dare = np.array([[Bd]])
Q_dare = np.array([[Qd]])
R_dare = np.array([[Rd]])

S = solve_discrete_are(A_dare, B_dare, Q_dare, R_dare)
S_scalar = S[0, 0] # Extract scalar value

# Calculate DLQR gain Kd
# Kd = (R + B'SB)^-1 B'SA
Kd_term1 = R_dare + B_dare.T @ S @ B_dare
Kd_term2 = B_dare.T @ S @ A_dare
Kd = np.linalg.inv(Kd_term1) @ Kd_term2
K_lqr = Kd[0, 0] # Extract scalar gain

# Check closed-loop stability for the linear system
Acl_lqr = Ad - Bd * K_lqr
print(f"DARE Solution S = {S_scalar:.3f} (for Q={Qd}, R={Rd})")
print(f"DLQR Gain K = {K_lqr:.3f}")
if abs(Acl_lqr) < 1:
    print(f"Linear Closed-Loop Pole = {Acl_lqr:.3f}. (STABLE)")
else:
    print(f"Linear Closed-Loop Pole = {Acl_lqr:.3f}. (UNSTABLE)")
# --- End DLQR Design ---

# Function to run the simulation (same as before, just uses K_lqr)
def run_simulation(x0_val, K_gain):
    t_history = [t_start]
    x_history = [x0_val]
    u_history = []
    current_t = t_start
    current_x = np.array([x0_val]) # State needs to be array

    for k in range(n_steps):
        x_deviation = current_x[0] - x_eq
        u_k = -K_gain * x_deviation
        u_history.append(u_k)

        def ode_interval(t, y):
            return np.array([nonlinear_system(y[0], u_k)])

        t_interval = (current_t, current_t + Ts)
        sol_interval = solve_ivp(
            ode_interval, t_interval, current_x, method='RK45',
            t_eval=[current_t + Ts]
        )

        if sol_interval.status != 0:
            print(f"Simulation failed at step {k} (t={current_t:.2f}) for x(0)={x0_val}")
            # Pad history...
            failed_steps = n_steps - k
            t_history.extend(np.linspace(current_t + Ts, t_end, failed_steps))
            x_history.extend([np.nan] * failed_steps)
            u_history.extend([np.nan] * failed_steps)
            break

        current_t += Ts
        current_x = sol_interval.y[:, -1]
        t_history.append(current_t)
        x_history.append(current_x[0])

        if abs(current_x[0]) > 10: # Stop if diverging
            print(f"State diverging at step {k} (t={current_t:.2f}) for x(0)={x0_val}. Stopping.")
            failed_steps = n_steps - k -1
            if failed_steps > 0:
                 t_history.extend(np.linspace(current_t + Ts, t_end, failed_steps))
                 x_history.extend([np.nan] * failed_steps)
            break

    return np.array(t_history), np.array(x_history), np.array(u_history)

# Run both simulations with DLQR gain
t_small, x_small, u_small = run_simulation(x0_small, K_lqr)
t_large, x_large, u_large = run_simulation(x0_large, K_lqr)

# 4. Plotting Results
plt.figure(figsize=(10, 8))
plt.subplot(2, 1, 1)
plt.plot(t_small, x_small, label=f'State x(t) (x0={x0_small})', linewidth=2)
plt.plot(t_large, x_large, label=f'State x(t) (x0={x0_large})', linewidth=2, linestyle='--')
plt.axhline(0, color='black', linewidth=0.5, linestyle=':')
plt.ylabel('State x')
plt.title(f'Scalar Nonlinear System with DLQR Controller (Q={Qd}, R={Rd})')
plt.grid(True)
plt.legend()
plt.ylim(min(np.nanmin(x_small), np.nanmin(x_large), -1) - 0.5, max(np.nanmax(x_small), np.nanmax(x_large), 1) + 0.5)

plt.subplot(2, 1, 2)
final_idx_small = len(u_small)
if final_idx_small > 0:
    # Use step plot for control signal visualization
    plt.step(t_small[:final_idx_small], u_small[:final_idx_small], where='post',
             label=f'Control u[k] (x0={x0_small})', linewidth=2)

final_idx_large = len(u_large)
if final_idx_large > 0:
    plt.step(t_large[:final_idx_large], u_large[:final_idx_large], where='post',
             label=f'Control u[k] (x0={x0_large})', linewidth=2, linestyle='--')

plt.ylabel('Control Input u')
plt.xlabel('Time (s)')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
