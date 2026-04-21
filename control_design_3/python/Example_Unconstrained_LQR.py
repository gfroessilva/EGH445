# Scenario 1: Unconstrained LQR for Nonzero Regulation (Damped 2nd Order System)

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_discrete_are, expm
from scipy.signal import cont2discrete
import matplotlib.pyplot as plt
import time

# print("--- Running Scenario 1: Unconstrained LQR (Nonzero Setpoint, Damped System) ---")

# --- 1. System Definition ---
# Continuous system: dx/dt = Ac*x + Bc*u
# x = [y, y_dot]^T = [x1, x2]^T
# d(x1)/dt = x2
# d(x2)/dt = x1 - c*x2 + u
c_damping = 0.5
Ac = np.array([[0., 1.], [1., -c_damping]])
Bc = np.array([[0.], [1.]])
n_states = Ac.shape[0]
m_inputs = Bc.shape[1]

def system_ode(t, x, u, Ac_mat, Bc_mat):
    u_val = u if np.isscalar(u) else u[0]
    dxdt = Ac_mat @ x + Bc_mat @ np.array([u_val])
    return dxdt.flatten()

# --- 2. Discretization ---
Ts = 0.1 # Sampling time
Ad, Bd, _, _, _ = cont2discrete((Ac, Bc, np.eye(n_states), np.zeros((n_states, m_inputs))), Ts, method='zoh')
# print("--- Discrete Model ---")
# print(f"Ad = \n{Ad}")
# print(f"Bd = \n{Bd}")

# --- 3. Setpoint ---
x_ref = np.array([1.0, 0.0]) # Target state y=1, y_dot=0
# Calculate corresponding u_ref for discrete system
# x_ref = Ad*x_ref + Bd*u_ref => Bd*u_ref = (I-Ad)*x_ref
# Use pseudo-inverse if Bd is not square: u_ref = pinv(Bd) @ (I-Ad) @ x_ref
try:
    u_ref_vec = np.linalg.pinv(Bd) @ (np.eye(n_states) - Ad) @ x_ref
    u_ref = u_ref_vec[0] # Extract scalar
    # Verify continuous u_ref = -x1_ref = -1
    # print(f"Setpoint: x_ref = {x_ref.T}, requires discrete u_ref = {u_ref:.3f} (cont. u_ref=-1.0)")
except Exception as e:
    print(f"Could not calculate discrete u_ref: {e}. Using continuous u_ref.")
    u_ref = -x_ref[0] # Use continuous calculation

# --- 4. DLQR Design (for error regulation) ---
# Weights for error cost: sum(xtilde' Q xtilde + utilde' R utilde)
Qd = np.diag([1.0, 0.1]) # Penalize position x1 more
Rd = np.array([[0.01]])   # Aggressive: Low penalty on control effort

try:
    S = solve_discrete_are(Ad, Bd, Qd, Rd)
    Kd_term1 = Rd + Bd.T @ S @ Bd
    Kd_term2 = Bd.T @ S @ Ad
    Kd = np.linalg.inv(Kd_term1) @ Kd_term2 # Gain is 1x2
    print(f"DLQR Gain for error K = {Kd}")
except Exception as e:
    print(f"DLQR Design failed: {e}"); exit()

# --- 5. Simulation Setup ---
t_start = 0
t_end = 6.0
n_steps = int(t_end / Ts)
x0 = np.array([0.0, 0.0]) # Initial condition (start from origin)

# --- Helper function for interval simulation ---
def simulate_interval(ode_func, t_start, x_start, duration, control_val):
    t_eval_interval = [t_start + duration]
    sol_interval = solve_ivp(
        ode_func, (t_start, t_start + duration), x_start, method='RK45',
        t_eval=t_eval_interval, args=(control_val,)
    )
    if sol_interval.status != 0: return None, False
    return sol_interval.y[:, -1], True

# --- 6. Simulation Loop ---
t_history = [t_start]
x_history = [x0]
u_history = []
current_t = t_start
current_x = x0.copy()
success_run = True

start_sim_time = time.time()
# print("Running simulation...")
for k in range(n_steps):
    x_error = current_x - x_ref
    # Control Law: u[k] = u_ref - K*(x[k]-x_ref)
    u_k_vec = u_ref - Kd @ x_error
    u_k = u_k_vec[0] # Extract scalar control
    # NO SATURATION APPLIED HERE
    u_history.append(u_k)

    def ode_lqr_interval(t, y, u_val):
        return system_ode(t, y, u_val, Ac, Bc)

    next_x, success = simulate_interval(ode_lqr_interval, current_t, current_x, Ts, u_k)

    if not success:
        print(f"ODE solver failed at t={current_t:.2f}")
        success_run = False; break

    current_t += Ts
    current_x = next_x
    t_history.append(current_t)
    x_history.append(current_x.copy())

    if np.max(np.abs(current_x)) > 100: # Divergence check
        print(f"State diverging excessively at t={current_t:.2f}. Stopping.")
        success_run = False; break

print(f"Simulation loop time: {time.time() - start_sim_time:.2f} s")

# --- 7. Plotting ---
if success_run:
    t_history = np.array(t_history)
    x_history = np.array(x_history).T # Shape is (n_states, n_points)
    u_history = np.array(u_history)
    n_plot = len(u_history)

    # Define consistent plot limits (adjust if needed)
    xlims = (t_start, t_end)
    ylims_x1 = (-0.5, 1.5)
    ylims_x2 = (-0.5, 1.5)
    ylims_u = (-4, 1) # Adjusted for u_ref=-1 and expected transient

    plt.figure("Unconstrained LQR (Nonzero Setpoint, Damped)", figsize=(10, 7))
    plt.rcParams.update({'font.size': 12})

    plt.subplot(3, 1, 1)
    plt.plot(t_history[:n_plot+1], x_history[0,:n_plot+1], label='$x_1$ (y)', lw=2)
    plt.axhline(x_ref[0], color='k', linestyle=':', label=f'Setpoint $x_1={x_ref[0]}$')
    plt.title(f'Unconstrained LQR (Setpoint x=[{x_ref[0]},{x_ref[1]}]^T)')
    plt.ylabel('Position $x_1$')
    plt.ylim(ylims_x1); plt.xlim(xlims); plt.grid(True); plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(t_history[:n_plot+1], x_history[1,:n_plot+1], label='$x_2$ (dy/dt)', lw=2)
    plt.axhline(x_ref[1], color='k', linestyle=':', label=f'Setpoint $x_2={x_ref[1]}$')
    plt.ylabel('Velocity $x_2$')
    plt.ylim(ylims_x2); plt.xlim(xlims); plt.grid(True); plt.legend()

    plt.subplot(3, 1, 3)
    plt.step(t_history[:n_plot], u_history[:n_plot], where='post', label='Control u[k]', lw=2)
    plt.axhline(u_ref, color='k', linestyle=':', label=f'Steady-State u={u_ref:.2f}')
    plt.ylabel('Control Input u')
    plt.xlabel('Time (s)')
    plt.ylim(ylims_u); plt.xlim(xlims); plt.grid(True); plt.legend()

    plt.tight_layout()
    plt.show()
else:
    print("Simulation failed, not plotting.")
