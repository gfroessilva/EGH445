% System Parameters
% \dot{x} = -ax + bu
a = 2;
b = 1;
x0 = 1; 
% Simulation Parameters
h = 0.8;  % Step size
t_end = 3.5; % End time
t = 0:h:t_end; % Time vector
t_real = 0:0.001:t_end; % Time for the analytical response

u = ones(size(t)); % Step input (u(t) = 1)

% Initialize solution vectors
x_fe = zeros(size(t));   % Forward Euler
x_be = zeros(size(t));   % Backward Euler
x_heun = zeros(size(t)); % Heun's Method
x_analytical = zeros(size(t)); % Analytical Solution

% Initial conditions
x_fe(1) = x0;
x_be(1) = x0;
x_heun(1) = x0;
x_analytical(1) = x0;

% Analytical Solution
for i = 1:length(t_real)-1
    x_analytical(i+1) = (x0 - 0.5)*exp(-a*t_real(i+1)) + 0.5;
end

% Numerical Integration
for i = 1:length(t)-1
    % Forward Euler
    x_fe(i+1) = (1 - h*a)*x_fe(i) + h*b*u(i);

    % Backward Euler
    x_be(i+1) = (x_be(i) + h*b*u(i+1)) / (1 + h*a);

    % Heun's Method
    x_tilde = x_heun(i) + h*(-a*x_heun(i) + b*u(i)); % Predictor
    x_heun(i+1) = x_heun(i) + (h/2)*((-a*x_heun(i) + b*u(i)) + (-a*x_tilde + b*u(i+1))); % Corrector

end

% Plotting
figure(1); clf
plot(t_real, x_analytical, 'k-', 'LineWidth', 4, 'DisplayName', 'Analytical'); hold on;
plot(t, x_fe, 'r--', 'LineWidth', 3, 'DisplayName', 'Forward Euler');
plot(t, x_be, 'b-.', 'LineWidth', 3, 'DisplayName', 'Backward Euler');
plot(t, x_heun, 'g:', 'LineWidth', 3, 'DisplayName', 'Heun''s Method');
hold off;

xlabel('Time (t)');
ylabel('x(t)');
title('Comparison of Numerical Integration Methods');
legend('Location', 'best');
grid on;

%% For a Second Order System Mass\-Spring\-Damper

% System Parameters (2nd Order System - Mass-Spring-Damper)
m = 1;   % Mass (kg)
k = 1;  % Spring constant (N/m) 
b = 40;   % Damping coefficient (Ns/m) 

% State-Space Matrices
A = [0 1; -k/m -b/m]
eigsA = eig(A);
disp(['The eigenvalues of the A matrix are [' num2str(eigsA') '].'])
disp(['The ratio between the real part of the eigenvalues is ' num2str(min(real(eigsA))/max(real(eigsA)))])
B = [0; 1/m];
C = [1 0];  % Output is position (you could change this)
D = 0;

% Initial Conditions
x0 = [1; 0];  % Initial position = 1, initial velocity = 0

%[text] Set the time step `h` and simulation time `t_end`:
% Simulation Parameters
h = 0.4;   % Step size
t_end = 15; % End time 
t = 0:h:t_end; % Time vector
u = zeros(size(t));  % Input:  Let's use a zero input for now (free response)
%u = ones(size(t)); % Step Input - uncomment to use a step input.
%u(t>=2) = 0;      % Pulse Input -  uncomment for a pulse input.

% Initialize solution vectors
x_fe = zeros(2, length(t));   % Forward Euler (2 states)
x_be = zeros(2, length(t));   % Backward Euler
x_heun = zeros(2, length(t)); % Heun's Method
x_analytical = zeros(2, length(t)); % Analytical Solution

% Initial conditions
x_fe(:, 1) = x0;
x_be(:, 1) = x0;
x_heun(:, 1) = x0;

% --- Numerical Integration ---
for i = 1:length(t)-1
    % Forward Euler
    x_fe(:, i+1) = (eye(2) + h*A)*x_fe(:, i) + h*B*u(i);

    % Backward Euler
    x_be(:, i+1) = (eye(2) - h*A) \ (x_be(:, i) + h*B*u(i+1)); % Note the matrix inverse

    % Heun's Method
    x_tilde = x_heun(:, i) + h*(A*x_heun(:, i) + B*u(i));    % Predictor
    x_heun(:, i+1) = x_heun(:, i) + (h/2)*((A*x_heun(:, i) + B*u(i)) + (A*x_tilde + B*u(i+1))); % Corrector
end

% --- Analytical Solution (Unified) ---
% 1. Characteristic Equation Roots
s = roots([m, b, k]);
s1 = s(1);
s2 = s(2);

% 2. Constants C1 and C2
C1 = (x0(2) - s2*x0(1)) / (s1 - s2);
C2 = (s1*x0(1) - x0(2)) / (s1 - s2);

% Check for NaN and handle critically damped case
if isnan(C1) || isnan(C2)
    C1 = x0(1);
    C2 = x0(2) + s1 * x0(1);
end

% 3. Calculate Analytical Solution
for i = 1:length(t)
    x_analytical(1, i) = real(C1*exp(s1*t(i)) + C2*exp(s2*t(i))); %Position
    x_analytical(2, i) = real(s1*C1*exp(s1*t(i)) + s2*C2*exp(s2*t(i))); %Velocity
end
% --- Plotting ---
figure(2);clf
subplot(2,1,1); % Plot position
plot(t, x_analytical(1,:), 'k-', 'LineWidth', 4, 'DisplayName', 'Analytical'); hold on;
plot(t, x_fe(1,:), 'r--', 'LineWidth', 3, 'DisplayName', 'Forward Euler');
plot(t, x_be(1,:), 'b-.', 'LineWidth', 3, 'DisplayName', 'Backward Euler');
plot(t, x_heun(1,:), 'g:', 'LineWidth', 3, 'DisplayName', 'Heun''s Method');
axis([0 t_end -1 1])
hold off;
xlabel('Time (t)');
ylabel('Position (x(t))');
title('Comparison of Numerical Integration Methods (Position)');
legend('Location', 'best');
grid on;

subplot(2,1,2); % Plot velocity
plot(t, x_analytical(2,:), 'k-', 'LineWidth', 4, 'DisplayName', 'Analytical'); hold on;
plot(t, x_fe(2,:), 'r--', 'LineWidth', 3, 'DisplayName', 'Forward Euler');
plot(t, x_be(2,:), 'b-.', 'LineWidth', 3, 'DisplayName', 'Backward Euler');
plot(t, x_heun(2,:), 'g:', 'LineWidth', 3, 'DisplayName', 'Heun''s Method');
axis([0 t_end -1 1])
hold off;
xlabel('Time (t)');
ylabel('Velocity (dx/dt)');
title('Comparison of Numerical Integration Methods (Velocity)');
legend('Location', 'best');
grid on;
