import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg
import control as ctrl
# Create the system matrices
G = np.array([[1, 0.1], [0, 0.8]])
H = np.array([[1, 0], [0, 1]])
C = np.array([[1, 0], [0, 1]]) # Outputs are the states
D = np.array([[0, 0], [0, 0]])

# Pole placement for MIMO system
K1 = np.array([[0.6, 0], [0, 0.2]]) 
K2 = np.array([[0.4, 0], [0, 0.4]]) 

# QLQR design for MIMO system
Q = np.array([[10, 0], [0, 10]]) # State weight matrix
R = np.array([[0.1, 0], [0, 0.1]]) # Control weight matrix

# Solve the Discrete Algebraic Riccati Equation (DARE)
# S is the solution to the DARE
S = linalg.solve_discrete_are(G, H, Q, R) 
# Calculate DLQR gain K
K3 = np.linalg.inv(R + H.T @ S @ H) @ (H.T @ S @ G)
print(f"K3 = [{K3[0, 0]:.3f}, {K3[0, 1]:.3f}; \n"
    f" {K3[1, 0]:.3f}, {K3[1, 1]:.3f}], "
    f" eig(G - H @ K3) = {np.linalg.eigvals(G - H @ K3)}")

# Compute the closed-loop system matrices
Gcl1 = G - H @ K1 
Gcl2 = G - H @ K2
Gcl3 = G - H @ K3

Ts = 0.1 # Sampling time (seconds)
# Create the closed-loop LTI systems
sys1 = ctrl.ss(Gcl1, H, C, D, Ts) 
sys2 = ctrl.ss(Gcl2, H, C, D, Ts) 
sys3 = ctrl.ss(Gcl3, H, C, D, Ts)

# Simulate for initial condition x0 = [1, 0]
x0 = np.array([1.5, 0.5])
t_end = 1.5 
t = np.arange(0, t_end, Ts) 
t1, y1 = ctrl.initial_response(sys1, T=t, X0=x0)
t2, y2 = ctrl.initial_response(sys2, T=t, X0=x0)
t3, y3 = ctrl.initial_response(sys3, T=t, X0=x0)

plt.figure(figsize=(10, 4))
plt.plot(t1, y1[0, :], label='K1 - x1', linewidth=3)
plt.plot(t1, y1[1, :], label='K1 - x2', linewidth=3)
plt.plot(t2, y2[0, :], label='K2 - x1', linestyle='--', 
    linewidth=3) 
plt.plot(t2, y2[1, :], label='K2 - x2', linestyle='--', 
    linewidth=3) 
plt.plot(t3, y3[0, :], label='K3 - x1', linestyle=':', 
    linewidth=3) 
plt.plot(t3, y3[1, :], label='K3 - x2', linestyle=':', 
    linewidth=3) 
plt.title('Initial Condition Response ($x_0=[1, 0]^T$) ' 
    'with Different K Matrices', fontsize=16) 
plt.xlabel('Time (steps * dt)', fontsize=14) 
plt.ylabel('State Values', fontsize=14) 
plt.legend(fontsize=12)
plt.grid(True) 
plt.show()
