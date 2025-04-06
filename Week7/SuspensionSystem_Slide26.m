%% Sinusoidal Tracking with state feedback controller and augmentation with internal model
clear
LineWidth = 5;
FontSize = 30;

%%% Parameters and System Definition
m = 1.5; k = 1; b = 1; y0 = 0;
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
fHz = 0.10; % frequency of the sinusoidal reference
omega = 2 * pi *fHz; % frequency of the sinusoidal reference
r(t >= 0) = .5 * sin(omega * t(t >= 0)); % sinusoidal reference

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
q1_history(:, 1) = 0; 
q2_history(:, 1) = 0;

% --- Define the internal model (oscillator) ---
Gq = [0 1; -1 2*cos(omega*Ts)]; % internal model of the reference
Hq = [0; 1]; % internal model of the reference

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the error state and the augmented system
e = [0 0]'; % initial error state

A_aug = [sysD_hat.A zeros(n_states, 2); -Hq*sysD_hat.C Gq];
B_aug = [sysD_hat.B; zeros(2, m_inputs)];
Br_aug = [zeros(n_states, 1); Hq];
C_aug = [sysD_hat.C zeros(1, 2)];
D_aug = sysD_hat.D;

sysD_aug = ss(A_aug, B_aug, C_aug, D_aug, Ts);

% Check if the augmented system is controllable
if rank(ctrb(sysD_aug.A, sysD_aug.B)) ~= size(sysD_aug.A, 1)
    error('The augmented system is not controllable!')
end
poles_desired = exp([-.5+1i -.5-1i -4 -3] * Ts); % Desired pole locations for the augmented system
K_aug = place(sysD_aug.A, sysD_aug.B, poles_desired); % Augmented state feedback gain
% K = place(sysD.A, sysD.B, [0.5 0.6]); % State feedback gain for the original system

% Calculate a feedforward control to follow the sinusoidal reference by inverting the DC gain of the system
Kf = (sysD.C*(eye(2) - (sysD.A-sysD.B*K_aug(1:n_states)))^(-1)*sysD.B)^(-1); 
u_ff = Kf * r; % (m x 1) % feedforward control to eliminate disturbance

% --- Simulation Loop ---
for k = 1:(N-1)
    % 1. Calculate Control Input u(k) based on state x(k)
    x_current = x_history(:, k);
    q1_current = q1_history(:, k);
    q2_current = q2_history(:, k);

    u_total = -K_aug(1:n_states) * x_current ...
              -K_aug(n_states+1:end) * [q1_current; q2_current];
    u_history(:, k) = u_total; % Store input u(k)

    % 2. Calculate Next State x(k+1)
    x_next = sysD.A * x_current + sysD.B * u_total + Bd_d * d(k);
    x_history(:, k+1) = x_next; % Store state x(k+1)
    
    % 3. Calculate Output y(k+1) (Optional)
    y_history(:, k+1) = C * x_next;
    q1_next = q2_current; % integral state update
    q2_next = -q1_current + 2*cos(omega*Ts)*q2_current + r(k) - C * x_current; 
    q1_history(:, k+1) = q1_next; % Store integral state q(k+1)
    q2_history(:, k+1) = q2_next; % Store integral state q(k+1)
end

f(10) = figure(10); f(10).Theme = 'light'; f(10).Color = 'white'; hold off
hold off
stairs(t, y_history, 'b', 'LineWidth', LineWidth);
hold on 
stairs(t, r, 'k--', 'LineWidth', LineWidth);
% stairs(t, d, 'r--', 'LineWidth', LineWidth);
legend('Wheel Position', 'Reference', 'Disturbance')
set(gca, 'FontSize', FontSize)
xlabel('Time (s)')
ylabel('Wheel Position (m)')
axis([0 t(end) -0.6 0.6])
grid on
hold off
