%% Sinusoidal Tracking with state feedback controller
clear
LineWidth = 5;
FontSize = 30;

%%% Parameters and System Definition
m = 1.5; k = 1; b = 1; y0 = 1;
A = [0 1; -k/m -b/m];
B = [0 1/m]';
C = [1 0];

sysC = ss(A,B,C,0);

Ts = 0.1;
sysD = c2d(sysC, Ts);

t = 0:Ts:10;
N = length(t);    % Number of simulation steps
x0 = [y0 0]';

Bd = [0 10/m]'; % how the disturbance affects the system

sysC_dist = ss(A, [B Bd], C, 0);
sysD_dist = c2d(sysC_dist, Ts);
Bd_d = sysD_dist.B(:, 2); % discrete-time disturbance input matrix

% Create a sinusoidal reference vector that starts at t=1 and ends at t_end
r = zeros(size(t));
fHz = 0.05; % frequency of the sinusoidal reference
omega = 2 * pi *fHz; % frequency of the sinusoidal reference
r(t >= 1) = .5 * sin(omega * t(t >= 1)); % sinusoidal reference

% Create a disturbance vector that starts at t=4 and ends at t=end
d = zeros(size(t));
d(t >= 4) = 0;

mhat = m; % estimated mass
khat = k; % estimated spring constant
bhat = b; % estimated damping coefficient
Ahat = [0 1; -khat/mhat -bhat/mhat];
Bhat = [0 1/mhat]';
sysC_hat = ss(Ahat, Bhat, C, 0);
sysD_hat = c2d(sysC_hat, Ts);

sysC_dist = ss(Ahat, [Bhat Bd], C, 0);
sysD_dist = c2d(sysC_dist, Ts);
Bd_hat_d = sysD_dist.B(:, 2); % discrete-time disturbance input matrix

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
q_history(:, 1) = 0; % initial error integral state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the error state and the augmented system
e = [0 0]'; % initial error state

A_aug = [sysD_hat.A zeros(n_states, p_outputs); -sysD_hat.C ones(p_outputs, p_outputs)];
B_aug = [sysD_hat.B; zeros(p_outputs, m_inputs)];
C_aug = [sysD_hat.C zeros(p_outputs, p_outputs)];
D_aug = sysD_hat.D;

sysD_aug = ss(A_aug, B_aug, C_aug, D_aug, Ts);

% Check if the augmented system is controllable
if rank(ctrb(sysD_aug.A, sysD_aug.B)) ~= size(sysD_aug.A, 1)
    error('The augmented system is not controllable!')
end
poles_desired = exp([-.5 -6 -15] * Ts); % Desired pole locations for the augmented system
K_aug = place(sysD_aug.A, sysD_aug.B, poles_desired); % Augmented state feedback gain
K = place(sysD.A, sysD.B, [0.5 0.6]); % State feedback gain for the original system

% Calculate a feedforward control to follow the sinusoidal reference by inverting the DC gain of the system
Kf = (sysD.C*(eye(2) - (sysD.A-sysD.B*K))^(-1)*sysD.B)^(-1); 
u_ff = Kf * r; % (m x 1) % feedforward control to eliminate disturbance

% --- Simulation Loop ---
for k = 1:(N-1)
    % 1. Calculate Control Input u(k) based on state x(k)
    x_current = x_history(:, k);
    q_current = q_history(:, k);

    u_total = -K(1:n_states) * x_current + u_ff(k);% - K_aug(n_states+1:end) * q_current;
    u_history(:, k) = u_total; % Store input u(k)

    % 2. Calculate Next State x(k+1)
    x_next = sysD.A * x_current + sysD.B * u_total + Bd_d * d(k);
    x_history(:, k+1) = x_next; % Store state x(k+1)
    
    % 3. Calculate Output y(k+1) (Optional)
    y_history(:, k+1) = C * x_next;
    q_next = q_current + r(k) - C * x_current; % integral state update
    q_history(:, k+1) = q_next; % Store integral state q(k+1)
end

f(8) = figure(8); f(8).Theme = 'light'; f(8).Color = 'white'; hold off
hold off
stairs(t, y_history, 'b', 'LineWidth', LineWidth);
hold on 
stairs(t, r, 'k--', 'LineWidth', LineWidth);
% stairs(t, d, 'r--', 'LineWidth', LineWidth);
legend('Wheel Position', 'Reference', 'Disturbance')
set(gca, 'FontSize', FontSize)
xlabel('Time (s)')
ylabel('Wheel Position (m)')
axis([0 t(end) -0.1 1.1])
grid on
hold off