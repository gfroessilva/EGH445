%% Zero state Regulation with disturbance
% The system is regulated using a state feedback controller with a non-zero initial state.
% The desired pole locations are defined, and the state feedback gain is calculated.
% The closed-loop system is then simulated with a step input.
% The simulation results are plotted to visualize the system response.

clear
LineWidth = 5;
FontSize = 30;

%%% Parameters and System Definition

m = 2; k = 1; b = 1; y0 = 1;
A = [0 1; -k/m -b/m];
B = [0 1/m]';
C = [1 0];

sysC = ss(A,B,C,0);

Ts = 0.1;
sysD = c2d(sysC, Ts);

t = 0:Ts:2.5;
N = length(t);    % Number of simulation steps
x0 = [y0 0]';

K = place(sysD.A, sysD.B, [0.5 0.6]);

Bd = [0 10/m]'; % how the disturbance affects the system
sysC_dist = ss(A, [B Bd], C, 0);

sysD_dist = c2d(sysC_dist, Ts);
Bd_d = sysD_dist.B(:, 2); % discrete-time disturbance input matrix
% Create a disturbance vector that starts at t=1 and ends at t=2.5
d = zeros(size(t));
d(t >= 1) = 1; % step disturbance at t=1
%%%%%%%%%%%%%%%%%%%%%


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
    x_next = sysD.A * x_current + sysD.B * u_total + Bd_d * d(k);
    x_history(:, k+1) = x_next; % Store state x(k+1)

    % 3. Calculate Output y(k+1) (Optional)
    y_history(:, k+1) = C * x_next;
end

f(2) = figure(2); f(2).Theme = 'light'; hold off
stairs(t, y_history, 'b', 'LineWidth', LineWidth);
hold on
stairs(t, d, 'k--', 'LineWidth', LineWidth);
legend('Wheel Position', 'Disturbance')
set(gca, 'FontSize', FontSize)
title('Zero-state Regulation with State Feedback Controller')
xlabel('Time (s)')
ylabel('Wheel Position (m)')
grid on
axis([0 t(end) -0.1 1.1])
hold off