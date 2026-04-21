import numpy as np
from scipy.integrate import solve_ivp
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

# Function to run the simulation for a given x0
def run_simulation_smooth(x0_val, points_per_interval=10):
    # Store results for smooth plotting
    t_history_smooth = [t_start]
    x_history_smooth = [x0_val]
    # Store results at discrete intervals for markers and control calc
    t_history_discrete = [t_start]
    x_history_discrete = [x0_val]
    u_history = [] # Control applied over the interval starting at t_history_discrete

    current_t = t_start
    current_x_start_of_interval = np.array([x0_val]) # State at beginning of interval

    # print(f"\nRunning smooth simulation for x(0) = {x0_val}...")
    for k in range(n_steps):
        # Calculate discrete control based on state at START of interval
        x_deviation = current_x_start_of_interval[0] - x_eq
        u_k = -K * x_deviation
        # Optional: Limit control effort
        # u_k = np.clip(u_k, -u_max, u_max)
        u_history.append(u_k)

        # Define ODE function for this interval (constant u_k)
        def ode_interval(t, y): # y is a 1-element array
            return np.array([nonlinear_system(y[0], u_k)])

        # Simulate one interval Ts, evaluating at multiple points
        t_eval_interval = np.linspace(current_t, current_t + Ts, points_per_interval, endpoint=True)
        sol_interval = solve_ivp(
            ode_interval,
            (current_t, current_t + Ts), # t_span for the interval
            current_x_start_of_interval, # Initial state for interval
            method='RK45',
            t_eval=t_eval_interval # Evaluate at these points
        )

        if sol_interval.status != 0:
            print(f"Simulation failed at step {k} (t={current_t:.2f}) for x(0)={x0_val}")
            # Pad history if needed
            # ... (padding logic can be added if needed) ...
            break

        # Update state for next interval START
        current_x_start_of_interval = sol_interval.y[:, -1] # State at the end
        current_t += Ts

        # Store history
        # Append points *after* the first one to avoid duplication
        t_history_smooth.extend(sol_interval.t[1:])
        x_history_smooth.extend(sol_interval.y[0, 1:])
        # Store discrete points
        t_history_discrete.append(current_t)
        x_history_discrete.append(current_x_start_of_interval[0])


        # Stop if state diverges excessively
        if abs(current_x_start_of_interval[0]) > 10:
            print(f"State diverging at step {k} (t={current_t:.2f}) for x(0)={x0_val}. Stopping.")
            break

    return (np.array(t_history_smooth), np.array(x_history_smooth), # Smooth results
            np.array(t_history_discrete), np.array(x_history_discrete), # Discrete results
            np.array(u_history)) # Control history

# --- Run both simulations ---
(t_small_smooth, x_small_smooth,
 t_small_discrete, x_small_discrete, u_small) = run_simulation_smooth(x0_small)

(t_large_smooth, x_large_smooth,
 t_large_discrete, x_large_discrete, u_large) = run_simulation_smooth(x0_large)


# --- MODIFIED Plotting Section ---
plt.figure(figsize=(10, 8))

# --- Plot states (Top Plot) ---
plt.subplot(2, 1, 1)
# Plot smooth trajectories
plt.plot(t_small_smooth, x_small_smooth,
         label=f'State x(t) (x0={x0_small})', linewidth=2)
plt.plot(t_large_smooth, x_large_smooth,
         label=f'State x(t) (x0={x0_large})', linewidth=2, linestyle='--')

plt.axhline(0, color='black', linewidth=0.5, linestyle=':')
plt.ylabel('State x')
plt.title('Scalar Nonlinear System with Linear Controller (Pole @ 0.5)')
plt.grid(True)
plt.legend()
# Adjust ylim automatically or set manually if needed
ylim_min = min(np.nanmin(x_small_smooth), np.nanmin(x_large_smooth), -1) - 0.5
ylim_max = max(np.nanmax(x_small_smooth), np.nanmax(x_large_smooth), 1) + 0.5
plt.ylim(ylim_min, ylim_max)

# --- Plot control inputs (Bottom Plot - using stairs) ---
plt.subplot(2, 1, 2)
# Ensure simulation ran successfully and generated history
final_idx_small = len(u_small)
if final_idx_small > 0:
    plt.stairs(u_small, t_small_discrete, label=f'Control u[k] (x0={x0_small})',
               baseline=None, linewidth=2, linestyle='--')

final_idx_large = len(u_large)
if final_idx_large > 0:
    plt.stairs(u_large[:final_idx_large-1], t_large_discrete, label=f'Control u[k] (x0={x0_large})',
               baseline=None, linewidth=2, linestyle='--')


plt.ylabel('Control Input u')
plt.xlabel('Time (s)')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
