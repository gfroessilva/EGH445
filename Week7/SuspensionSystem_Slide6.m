% This script simulates a simple suspension system using state-space representation.
% It converts a continuous-time state-space model to a discrete-time model.
% The system parameters are defined, and the state-space matrices are constructed.
% The continuous-time system is then converted to a discrete-time system with a specified sampling time.

clear
%% Plotting configuration
LineWidth = 5;
FontSize = 30;

%% Parameters and System Definition
m = 2; k = 1; b = 1; y0 = 1;
A = [0 1; -k/m -b/m];
B = [0 1/m]';
C = [1 0];

sysC = ss(A,B,C,0);

Ts = 0.1;
sysD = c2d(sysC, Ts);

%% Zero-state Regulation
% The system is regulated using a state feedback controller.
% The desired pole locations are defined, and the state feedback gain is calculated.
% The closed-loop system is then simulated with a step input.
% The simulation results are plotted to visualize the system response.

K = place(sysD.A, sysD.B, [0.5 0.6]);

t = 0:Ts:2.5;
N = length(t);    % Number of simulation steps
x0 = [y0 0]';
% sysD_CL = ss(sysD.A - sysD.B*K, sysD.B, sysD.C, sysD.D, Ts);

% --- Initialize History Arrays ---
n_states = size(sysD.A, 1);
m_inputs = size(sysD.B, 2);
p_outputs = size(sysD.C, 1);

x_history = zeros(n_states, N);
u_history = zeros(m_inputs, N-1); % Control input from k=0 to N-2
y_history = zeros(p_outputs, N);
% --- Set Initial Condition ---
x_history(:, 1) = x0;
y_history(:, 1) = sysD.C * x0;

% --- Simulation Loop ---
for k = 1:(N-1)
    % 1. Calculate Control Input u(k) based on state x(k)
    x_current = x_history(:, k);
    u_total = -K * x_current;
    u_history(:, k) = u_total; % Store input u(k)

    % 2. Calculate Next State x(k+1)
    x_next = sysD.A * x_current + sysD.B * u_total;
    x_history(:, k+1) = x_next; % Store state x(k+1)

    % 3. Calculate Output y(k+1) (Optional)
    y_history(:, k+1) = C * x_next;
end

f(1) = figure(1); f(1).Theme = 'light'; hold off
stairs(t, y_history, 'b', 'LineWidth', LineWidth);
legend('Wheel Position')
set(gca, 'FontSize', FontSize)
title('Zero-state Regulation with State Feedback Controller')
xlabel('Time (s)')
ylabel('Wheel Position (m)')
grid on
axis([0 t(end) -0.1 1.1])