% System Matrix (Stable)
A = [-2 1; -1 -3];

% --- Eigenvalue Check (for comparison) ---
eigenvalues = eig(A);
disp('Eigenvalues:');
disp(eigenvalues);
if all(real(eigenvalues) < 0)
    disp('System is asymptotically stable (Eigenvalues).');
end

% --- Lyapunov Equation ---
Q = eye(2);  % Choose Q = I (positive definite)
P = lyap(A, Q);
disp('Eigenvalues of P:');
disp(eig(P));
if all(eig(P) > 0)
    disp('System is asymptotically stable (Lyapunov).');
end

% --- Simulation ---
x0 = [2; -1];   % Initial condition
t_end = 5;      % Simulation time
h = 0.01;       % Time step
t = 0:h:t_end;
x = zeros(2, length(t));
x(:, 1) = x0;

% Use the closed-form solution (matrix exponential) for accuracy
for i = 1:length(t) - 1
    x(:, i+1) = expm(A*h) * x(:, i);
end

% --- Lyapunov Function Values ---
V = zeros(1, length(t));
for i = 1:length(t)
    V(i) = x(:, i)' * P * x(:, i);
end

% --- Create a Meshgrid for Contour Plot ---
[X1, X2] = meshgrid(-3:0.1:3, -3:0.1:3);
V_grid = zeros(size(X1));
for i = 1:size(X1, 1)
    for j = 1:size(X1, 2)
        x_grid = [X1(i, j); X2(i, j)];
        V_grid(i, j) = x_grid' * P * x_grid;
    end
end

% --- Create the Animation ---
figure;

% --- Subplot 1: State Trajectory and Lyapunov Contours ---
subplot(2, 1, 1);
contour(X1, X2, V_grid, 20); % Plot contours of V(x)
hold on;
h_traj = plot(x(1, 1), x(2, 1), 'r-', 'LineWidth', 2); % Trajectory
h_point = plot(x(1, 1), x(2, 1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Current state
hold off;
xlabel('x1');
ylabel('x2');
title('State Trajectory and Lyapunov Function Contours');
axis equal;  % Important for correct visualization of contours
grid on;

% --- Subplot 2: Lyapunov Function Value over Time ---
subplot(2, 1, 2);
h_V = plot(t(1), V(1), 'b-', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('V(x(t))');
title('Lyapunov Function Value');
grid on;
ylim([0, max(V)*1.1]); % Set y-axis limits

% --- Animation Loop ---
for i = 1:length(t)
    % Update state trajectory plot
    set(h_traj, 'XData', x(1, 1:i), 'YData', x(2, 1:i));
    set(h_point, 'XData', x(1, i), 'YData', x(2, i));

    % Update Lyapunov function plot
    set(h_V, 'XData', t(1:i), 'YData', V(1:i));

    drawnow;  % Force MATLAB to update the plot
    pause(0.01); % Control animation speed
    
     % Capture the frame
    frame = getframe(gcf); %Get the current frame
    im = frame2im(frame);%transform the frame to image

end