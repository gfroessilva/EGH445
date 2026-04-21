import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_discrete_are

# Iterate over different Q values and simulate the system
Q_values = [1, 10, 50, 100]  # Different state penalty values
R = 0.1  # Fixed control penalty
Ad = 0.905
Bd = 0.095
Ts = 0.1
x0 = 2  # Initial condition

# Store results for plotting
results = []

for Q in Q_values:
  # Solve DARE for each Q
  S = solve_discrete_are(np.array([[Ad]]), np.array([[Bd]]), np.array([[Q]]), np.array([[R]]))
  K = np.linalg.inv(R + Bd**2 * S) @ (Bd * S * Ad)
  K = K[0, 0]  # Extract scalar gain

  # Simulate the system
  t_history = [0]
  x_history = [x0]
  current_x = x0
  for k in range(50):  # Simulate for 50 steps
    u = -K * current_x
    current_x = Ad * current_x + Bd * u
    t_history.append((k + 1) * Ts)
    x_history.append(current_x)

  results.append((Q, t_history, x_history))

# Plot results
plt.figure(figsize=(10, 6))
for Q, t_history, x_history in results:
  plt.plot(t_history, x_history, label=f"Q={Q}", linewidth=2)

plt.title("DLQR State Response for Different Q Values")
plt.xlabel("Time (s)")
plt.ylabel("State x")
plt.legend()
plt.grid(True)
plt.show()
