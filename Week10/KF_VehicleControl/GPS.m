% This example shows a vehicle being controlled towards waypoints using
% only GPS measurements. The states are estimated by a Kalmna Filter
%
% Adapted from: Prof Aurelio Tergolina Salton
% PUCRS, Brazil, 2015 
% Last modified: Dr Guilherme Froes Silva, May 2025

clear all;
clc;

%% Simulation Parameters
samplingPeriod = 0.5; % Sampling period (T)
% N_max_steps = 50; % Example: Intended maximum number of simulation steps (currently unused as loop duration is determined by waypoints)

% Initialize figure for map display
figure(1); clf; % Clear figure
mapImage = imread('mapa.png'); % Load map image
imshow(mapImage); % Display map
hold on; % Keep map for subsequent plots

% Noise levels
gpsNoiseLevel = 100;       % Noise magnitude for GPS position measurements
accelerometerNoiseLevel = 0.001; % Noise magnitude for accelerometer measurements

% Waypoint navigation parameters
waypointProximityThresholds = [10 10 10 10 10 10 10 300 1000]; % Proximity threshold to consider a waypoint reached
targetWaypoints = [330 340 400 405 499 513 514  30 -1000   0;  % x-coordinates of waypoints
                   450 410 405 435 447 320 274 275  1000   0]; % y-coordinates of waypoints

%% System Definition (State-Space Model)
% Continuous-time state-space matrices
% State vector: [x_position; x_velocity; y_position; y_velocity]
A_c = [0 1 0 0;  % Continuous-time state matrix
       0 0 0 0;
       0 0 0 1;
       0 0 0 0];
% Discretized state transition matrix (A_d or Phi_k)
% Since A_c^2 = 0 for this system, expm(A_c*T) = I + A_c*T.
% The formula used, I + A_c*T + A_c^2*T^2/2, is exact here.
A_d = eye(4) + samplingPeriod*A_c + (A_c^2 * samplingPeriod^2)/2;

B_c = [0 0;  % Continuous-time input matrix
       1 0;  % Assuming input u = [ax; ay] (accelerations)
       0 0;
       0 1];
% Discretized input matrix (B_d or Gamma_k)
% B_d = samplingPeriod*B_c is a common approximation (Euler method for integral of B).
% Exact for this system would be (I*samplingPeriod + A_c*samplingPeriod^2/2)*B_c.
B_d = samplingPeriod*B_c;

C_obs = [1 0 0 0;  % Observation matrix (measures x_pos, y_pos)
         0 0 1 0];
C_d = C_obs; % Discrete observation matrix (assumed same as continuous)

% Kalman Filter Covariance Matrices
% Q_kf: Measurement noise covariance matrix (R_k in some standard KF notations)
% Represents uncertainty in sensor measurements (GPS). Diagonal 2x2 matrix.
Q_cov_measurement = eye(2) * gpsNoiseLevel^2; % Assuming gpsNoiseLevel is std dev, if magnitude use as is or adjust. Original was eye(2)*50.
% R_kf: Process noise covariance matrix (Q_k in some standard KF notations)
% Represents uncertainty in the dynamic model. Diagonal 4x4 matrix.
R_cov_process = eye(4) * accelerometerNoiseLevel^2; % Assuming accelerometerNoiseLevel is std dev. Original was eye(4)*0.001.

% P_initial: Initial state estimation error covariance
% Represents uncertainty in the initial state estimate. Diagonal 4x4 matrix.
P_est_cov = eye(4) * 500000;

%% Initial Conditions
x_estimated = [0 0 0 0]';       % Initial state estimate for the Kalman filter [x_pos, x_vel, y_pos, y_vel]'
x_true = [300 0 500 0]';        % Initial true state of the system [x_pos, x_vel, y_pos, y_vel]'
u_control = [0 0]';             % Initial control input [ax_command, ay_command]'

% Plot waypoints on the map
plot(targetWaypoints(1,:), targetWaypoints(2,:), 'ko', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'y');
% legend('Waypoints');

%% Simulation Loop
currentWaypointIdx = 1; % Index for the current target waypoint
simStep = 0;            % Simulation step counter

% Preallocate arrays for storing history (optional, for efficiency)
max_duration_steps = 100000; % Or estimate based on waypoints
x_true_history = zeros(4, max_duration_steps+1);
x_estimated_history = zeros(4, max_duration_steps+1);
z_measurements_history = zeros(2, max_duration_steps);
u_control_history = zeros(2, max_duration_steps+1);
x_true_history(:,1) = x_true;
x_estimated_history(:,1) = x_estimated;
u_control_history(:,1) = u_control;


% Loop while there are waypoints to navigate to
% Note: Original loop `j<length(r)` processes waypoints 1 to length(r)-1.
% If length(r) is 10, it processes waypoints 1 through 9.
while currentWaypointIdx < size(targetWaypoints, 2)
    simStep = simStep + 1;

    % Example: Dynamically change noise covariance (currently commented out)
    % if (simStep == 200)
    %     R_cov_measurement = eye(2)*30000; % Increase measurement noise
    % end

    % Check if the vehicle has reached the current waypoint
    % Uses true position at the beginning of the current step (x_true(:,simStep))
    current_pos_true = [x_true(1,simStep); x_true(3,simStep)];
    target_pos = [targetWaypoints(1,currentWaypointIdx); targetWaypoints(2,currentWaypointIdx)];
    distance_to_target = norm(target_pos - current_pos_true, 2);

    if distance_to_target < waypointProximityThresholds(currentWaypointIdx)
        disp(['Reached waypoint: ', num2str(currentWaypointIdx)]);
        % keyboard; % Pause for debugging if waypoint is reached
        currentWaypointIdx = currentWaypointIdx + 1;
        if currentWaypointIdx > size(targetWaypoints, 2) % All waypoints processed
             disp('All waypoints reached.');
             break; 
        end
    end

    % Simulate the true system dynamics (e.g., using XYmass function)
    % x_true(:, simStep+1) = XYmass(x_true(:, simStep), samplingPeriod, u_control(:, simStep));
    % Using state-space model directly if XYmass is not defined:
    x_true(:, simStep+1) = A_d * x_true(:, simStep) + B_d * u_control(:, simStep); % Process noise could be added here for more realism

    % Generate sensor measurements (GPS)
    % Measurement z = C*x_true + noise
    % Noise is uniform in [-0.5, 0.5] * gpsNoiseLevel
    measurement_noise = gpsNoiseLevel * (rand(2,1) - 0.5);
    z_measurements(:, simStep) = C_d * x_true(:, simStep+1) + measurement_noise;  % [x_pos_measured; y_pos_measured]

    % Generate noisy control input measurement (e.g., from accelerometer)
    % Noise is uniform in [-0.5, 0.5] * accelerometerNoiseLevel, applied independently to each axis
    control_input_noise = accelerometerNoiseLevel * (rand(2,1) - 0.5);
    u_measured(:, simStep) = u_control(:, simStep) + control_input_noise;

    % Kalman Filter Step
    % [x_est_next, P_cov_next] = Fkalman(x_est_current, P_cov_current, A_d, B_d, C_d, T, Q_process, R_measurement, z_current, u_measured_current)
    % Ensure Q_cov_process and R_cov_measurement are correctly passed based on Fkalman's expected signature.
    % Original call: Fkalman(x(:,k),P,A_d,B_d,C_d,T,R,Q,z(:,k),w(:,k))
    % R (original) is process noise Q_cov_process
    % Q (original) is measurement noise R_cov_measurement
    [x_estimated(:, simStep+1), P_est_cov] = Fkalman(x_estimated(:, simStep), P_est_cov, A_d, B_d, C_d, samplingPeriod, R_cov_process, Q_cov_measurement, z_measurements(:, simStep), u_measured(:, simStep));

    % Control Law (e.g., using 'controle' function)
    % u_control_next = computeControlInput(state_for_control, target_waypoint)
    % Original used: controle(x_estimated_next, targetWaypoint)
    if currentWaypointIdx <= size(targetWaypoints, 2) % Ensure waypoint index is valid
        % u_control(:, simStep+1) = controle(x_estimated(:, simStep+1), targetWaypoints(:, currentWaypointIdx));
        u_control(:, simStep+1) = controle(x_estimated(:, simStep+1), targetWaypoints(:, currentWaypointIdx));
        % u_control(:, simStep+1) = controle(x_true(:, simStep+1), targetWaypoints(:, currentWaypointIdx));
    else
        u_control(:, simStep+1) = [0;0]; % No target, stop
    end
    
    % Live Plotting
    if simStep == 1 % First step, plot initial true and estimated positions
        h_true_path = plot(x_true(1,1:simStep+1), x_true(3,1:simStep+1), 'k-', 'LineWidth', 2); % Plot initial true position and path start
        h_gps_meas = plot(z_measurements(1,simStep), z_measurements(2,simStep), 'r.', 'MarkerSize', 15); % Current GPS measurement
        
        h_est_pos = plot(x_estimated(1,simStep+1), x_estimated(3,simStep+1), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b'); % Current estimated position        
    elseif simStep == 2
        h_est_path = plot(x_estimated(1,1:simStep+1), x_estimated(3,1:simStep+1), 'b-', 'LineWidth', 2); % Plot initial estimated position and path start
    else
        % Update path data
        set(h_true_path, 'XData', x_true(1,1:simStep+1), 'YData', x_true(3,1:simStep+1));
        set(h_gps_meas, 'XData', z_measurements(1,simStep), 'YData', z_measurements(2,simStep));
        set(h_est_path, 'XData', x_estimated(1,3:simStep+1), 'YData', x_estimated(3,3:simStep+1));
        set(h_est_pos, 'XData', x_estimated(1,simStep+1), 'YData', x_estimated(3,simStep+1));
    end
    
    

    title(['Simulation Step: ', num2str(simStep)]);
    drawnow;
    pause(0.02); % Optional pause for slower animation
end
disp('Simulation finished.');

%% Plot Results
finalSimStep = simStep; % Number of steps actually run

% Time vectors for plotting
time_measurements = (1:finalSimStep) * samplingPeriod;
time_states = (0:finalSimStep) * samplingPeriod; % Includes initial state at t=0


figure(2); clf;

% Subplot 1: Positions (True, Measured, Estimated)
subplot(2,1,1);
hold on;
plot(time_measurements, z_measurements(1,1:finalSimStep), 'r.', 'DisplayName', 'x Measured (GPS)');
plot(time_measurements, z_measurements(2,1:finalSimStep), 'm.', 'DisplayName', 'y Measured (GPS)');
plot(time_states, x_true(1,1:finalSimStep+1), 'k-', 'LineWidth', 1.5, 'DisplayName', 'x True');
plot(time_states, x_true(3,1:finalSimStep+1), 'color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'DisplayName', 'y True');
plot(time_states, x_estimated(1,1:finalSimStep+1), 'b--', 'LineWidth', 1.5, 'DisplayName', 'x Estimated');
plot(time_states, x_estimated(3,1:finalSimStep+1), 'c--', 'LineWidth', 1.5, 'DisplayName', 'y Estimated');
hold off;
title('True, Measured, and Estimated Positions');
xlabel('Time (s)');
ylabel('Position (units)');
legend('show', 'Location', 'best');
grid on;

% Subplot 2: Velocities (True, Estimated)
subplot(2,1,2);
hold on;
plot(time_states, x_true(2,1:finalSimStep+1), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Vx True');
plot(time_states, x_true(4,1:finalSimStep+1), 'color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'DisplayName', 'Vy True');
plot(time_states, x_estimated(2,1:finalSimStep+1), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Vx Estimated');
plot(time_states, x_estimated(4,1:finalSimStep+1), 'c--', 'LineWidth', 1.5, 'DisplayName', 'Vy Estimated');
hold off;
title('True and Estimated Velocities');
xlabel('Time (s)');
ylabel('Velocity (units/s)');
legend('show', 'Location', 'best');
grid on;
