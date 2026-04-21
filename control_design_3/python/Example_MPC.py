# Scenario 3: MPC for Nonzero Regulation (Damped 2nd Order System)

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_discrete_are, expm
from scipy.signal import cont2discrete
import cvxpy as cp # Requires cvxpy and a solver like OSQP
import matplotlib.pyplot as plt
import time

# print("--- Running Scenario 3: MPC (Nonzero Setpoint, Damped System) ---")

# --- 1. System Definition ---
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
Ts = 0.1
Ad, Bd, _, _, _ = cont2discrete((Ac, Bc, np.eye(n_states), np.zeros((n_states, m_inputs))), Ts, method='zoh')
# print("--- Discrete Model ---")
# print(f"Ad = \n{Ad}")
# print(f"Bd = \n{Bd}")

# --- 3. Setpoint ---
x_ref = np.array([[1.0], [0.0]]) # Target state y=1, y_dot=0 (use column vector shape)
try:
    u_ref_vec = np.linalg.pinv(Bd) @ (np.eye(n_states) - Ad) @ x_ref
    u_ref = u_ref_vec[0,0] # Extract scalar
    # print(f"Setpoint: x_ref = {x_ref.T}, requires discrete u_ref = {u_ref:.3f} (cont. u_ref=-1.0)")
except Exception as e:
    print(f"Could not calculate discrete u_ref: {e}. Using continuous u_ref.")
    u_ref = -x_ref[0,0] # Use continuous calculation

# --- 4. MPC Setup ---
P = 30 # Prediction Horizon
# Use same base weights as LQR example, but consider tuning R_mpc
Q_mpc = np.diag([1.0, 0.1])
# R_mpc = np.array([[0.01]]) # Original aggressive LQR R
R_mpc = np.array([[0.1]]) # Use less aggressive R for MPC? Try this.
# Calculate terminal weight Qf from DARE solution S using MPC weights
try:
    S = solve_discrete_are(Ad, Bd, Q_mpc, R_mpc)
    Qf_mpc = S
    # print(f"Using terminal weight Qf=S based on MPC Q/R: \n{Qf_mpc}")
except Exception as e:
    print(f"DARE solve failed ({e}), using Qf=Q"); Qf_mpc = Q_mpc

# Constraints
u_min = -2.5
u_max = 0.5
print(f"Control Limits: [{u_min}, {u_max}]")

# --- 5. CVXPY Optimization Problem Setup ---
U = cp.Variable((m_inputs, P), name='U')
x_k_param = cp.Parameter(n_states, name='x_k_param') # Current state parameter
cost = 0.0
constraints = []
x_pred = x_k_param # Use parameter for current state

for t in range(P):
    # Penalize deviation from reference state & absolute control input
    cost += cp.quad_form(x_pred - x_ref.flatten(), Q_mpc) + cp.quad_form(U[:, t], R_mpc)
    x_next = Ad @ x_pred + Bd @ U[:, t]
    constraints += [U[:, t] >= u_min, U[:, t] <= u_max]
    x_pred = x_next
# Terminal cost penalizes deviation from x_ref
cost += cp.quad_form(x_pred - x_ref.flatten(), Qf_mpc)

objective = cp.Minimize(cost)
problem = cp.Problem(objective, constraints)
print("MPC CVXPY Problem defined.")

# --- 6. Simulation Setup ---
t_start = 0
t_end = 6.0
n_steps = int(t_end / Ts)
x0 = np.array([0.0, 0.0]) # Initial condition

# --- Helper function for interval simulation ---
def simulate_interval(ode_func, t_start, x_start, duration, control_val):
    t_eval_interval = [t_start + duration]
    sol_interval = solve_ivp(
        ode_func, (t_start, t_start + duration), x_start, method='RK45',
        t_eval=t_eval_interval, args=(control_val,)
    )
    if sol_interval.status != 0: return None, False
    return sol_interval.y[:, -1], True

# --- 7. Simulation Loop ---
t_history = [t_start]
x_history = [x0]
u_history = [] # Applied MPC control
current_t = t_start
current_x = x0.copy()
success_run = True
infeasibility_count = 0

start_sim_time = time.time()
# print("Running simulation...")
for k in range(n_steps):
    x_k_param.value = current_x # Update parameter
    u_k = u_ref # Default if solver fails
    try:
        problem.solve(solver=cp.OSQP, warm_start=True, verbose=False)
        if problem.status == cp.OPTIMAL or problem.status == cp.OPTIMAL_INACCURATE:
            u_k = U.value[0, 0]
        else:
            print(f"MPC Warning: Solver status {problem.status} at k={k}")
            infeasibility_count += 1
            if U.value is not None: u_k = U.value[0, 0]
            if problem.status == cp.INFEASIBLE:
                 print("MPC Problem is infeasible. Stopping."); success_run = False; break
    except Exception as e:
        print(f"MPC Solver failed at k={k}: {e}"); infeasibility_count += 1

    # Apply saturation as safety/clip numerical tolerances
    u_k_sat = np.clip(u_k, u_min, u_max)
    u_history.append(u_k_sat)

    def ode_mpc_interval(t, y, u_val):
        return system_ode(t, y, u_val, Ac, Bc)

    next_x, success = simulate_interval(ode_mpc_interval, current_t, current_x, Ts, u_k_sat)

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

print(f"MPC Simulation loop time: {time.time() - start_sim_time:.2f} s")
if infeasibility_count > 0: print(f"Warning: MPC Solver non-optimal status {infeasibility_count} times.")

# --- 8. Plotting ---
if success_run:
    t_history = np.array(t_history)
    x_history = np.array(x_history).T
    u_history = np.array(u_history)
    n_plot = len(u_history)

    # Use consistent plot limits from Block 1 & 2 (or adjust)
    xlims = (t_start, t_end)
    ylims_x1 = (-0.5, 1.5)
    ylims_x2 = (-0.5, 1.5)
    ylims_u = (-4, 1) # Use same limits

    plt.figure("MPC (Nonzero Setpoint, Damped)", figsize=(10, 7))
    plt.rcParams.update({'font.size': 12})

    plt.subplot(3, 1, 1)
    plt.plot(t_history[:n_plot+1], x_history[0,:n_plot+1], label='$x_1$ (y)', lw=2)
    plt.axhline(x_ref[0], color='k', linestyle=':', label=f'Setpoint $x_1={x_ref[0,0]}$')
    plt.title(f'MPC Response (Setpoint x=[{x_ref[0,0]},{x_ref[1,0]}]^T, u limits=[{u_min},{u_max}], P={P})')
    plt.ylabel('Position $x_1$')
    plt.ylim(ylims_x1); plt.xlim(xlims); plt.grid(True); plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(t_history[:n_plot+1], x_history[1,:n_plot+1], label='$x_2$ (dy/dt)', lw=2)
    plt.axhline(x_ref[1], color='k', linestyle=':', label=f'Setpoint $x_2={x_ref[1,0]}$')
    plt.ylabel('Velocity $x_2$')
    plt.ylim(ylims_x2); plt.xlim(xlims); plt.grid(True); plt.legend()

    plt.subplot(3, 1, 3)
    plt.step(t_history[:n_plot], u_history[:n_plot], where='post', label='Applied MPC u[k]', lw=2)
    plt.axhline(u_ref, color='k', linestyle=':', label=f'Steady-State u={u_ref:.2f}')
    plt.axhline(u_max, color='r', linestyle=':', label='Limit')
    plt.axhline(u_min, color='r', linestyle=':')
    plt.ylabel('Control Input u')
    plt.xlabel('Time (s)')
    plt.ylim(ylims_u); plt.xlim(xlims); plt.grid(True); plt.legend()

    plt.tight_layout()
    plt.show()
else:
    print("Simulation failed, not plotting.")
