% Main Simulation Script for Vehicle Navigation and Kalman Filtering
% Purpose: To simulate a 2D vehicle navigating waypoints, generate noisy
%          sensor readings (GPS for position, accelerometer for input),
%          and use a Kalman filter to estimate the vehicle's true state
%          (position and velocity).
% Author: Guilherme Froes Silva (adapted for classroom use)
% Date: 07 May 2025

clear all; % Clear all variables from workspace
clc;       % Clear the command window

%% Simulation Parameters
disp('Setting up simulation parameters...');
samplingTime = 0.5; % Simulation time step (seconds)
% N_steps = 50;    % Total number of simulation steps (original, but loop is waypoint-driven)

% --- Figure for Real-time Animation ---
figure(1); clf; % Create or clear figure 1
try
    mapImage = imread('mapa.png'); % Load a background map image
    imshow(mapImage);              % Display the map
    hold on;                       % Hold the plot for overlaying trajectories
    disp('Map loaded successfully.');
catch
    disp('Map file "mapa.png" not found. Plotting on a blank background.');
    axis([0 1000 0 1000]); % Define plot limits if map is not available
    hold on;
end

% --- Noise Parameters ---
gpsNoiseMagnitude = 50;    % Standard deviation of GPS position noise (units)
accelNoiseMagnitude = 0.001; % Standard deviation of accelerometer noise (units/s^2)

% --- Waypoint Navigation ---
% 'waypoints' stores the [x; y] coordinates of target locations
waypoints = [330 340 400 405 499 513 514   30 -1000   0;  % x-coordinates
             450 410 405 435 447 320 274  275  1000   0]; % y-coordinates
% 'waypointProximityThreshold' defines how close the vehicle needs to be to a waypoint to consider it "reached"
waypointProximityThreshold = [5 5 1 6 5 1 5 330 1000];

%% System Model (for Kalman Filter design)
% State vector: x_state = [x_position; x_velocity; y_position; y_velocity]
disp('Defining system model matrices...');
% Continuous-time state matrix (double integrator model)
A_continuous = [0 1 0 0;
                0 0 0 0;
                0 0 0 1;
                0 0 0 0];

% Continuous-time input matrix (maps acceleration input to velocity change)
B_continuous = [0 0;
                1 0;
                0 0;
                0 1];

% Discretisation of system matrices (using second-order Taylor expansion for A)
% x_k+1 = A_discrete * x_k + B_discrete * u_k
A_discrete = eye(4) + samplingTime * A_continuous + (samplingTime^2 / 2) * (A_continuous^2);
B_discrete = samplingTime * B_continuous;

% Output matrix (measures x_position and y_position)
% z_measurement = C_discrete * x_state
C_discrete = [1 0 0 0;
              0 0 1 0];

%% Kalman Filter Covariance Matrices
disp('Setting up Kalman Filter covariance matrices...');
% Measurement Noise Covariance (R_kf in many texts, Q here)
% Represents the uncertainty in the sensor measurements (GPS).
% Assumes uncorrelated noise between x and y position readings.
Q_measurementNoiseCov = eye(2) * gpsNoiseMagnitude^2; % Variance = sigma^2

% Process Noise Covariance (Q_kf in many texts, R here)
% Represents the uncertainty in the system model (e.g., unmodelled forces).
R_processNoiseCov = eye(4) * 0.001; % Small values indicate high confidence in the model

% Initial Estimate Error Covariance
% Represents the initial uncertainty in the filter's state estimate.
P_estimateErrorCov = eye(4) * 50000; % Large values indicate high initial uncertainty

%% Initial Conditions
disp('Setting initial conditions...');
% Filter's initial state estimate [px, vx, py, vy]'
initialStateEstimate = [0 0 0 0]';

% System's true initial state [px, vx, py, vy]'
trueInitialState = [300 0 500 0]';

% Initial control input (accelerations [ax; ay]')
currentControlInput = [0 0]';

% Pre-allocate storage for simulation data (optional, but good practice for speed)
% max_steps = length(waypointProximityThreshold) * 100; % Estimate max steps
% trueStates = zeros(4, max_steps);
% estimatedStates = zeros(4, max_steps);
% measurements = zeros(2, max_steps);
% controlInputs = zeros(2, max_steps);
% trueStates(:,1) = trueInitialState;
% estimatedStates(:,1) = initialStateEstimate;

%% Main Simulation Loop
disp('Starting simulation loop...');
currentWaypointIndex = 1;
timeStep = 0;

% Pointers to current states/inputs for cleaner loop code
estimatedState = initialStateEstimate;
trueState = trueInitialState;

while currentWaypointIndex <= length(waypointProximityThreshold)
    timeStep = timeStep + 1;

    % --- Dynamic Event: Simulate sudden increase in sensor noise ---
    % if timeStep == 200
    %     disp('Simulating increased GPS noise at timestep 200.');
    %     Q_measurementNoiseCov = eye(2) * 30000; % Drastically increase measurement noise
    % end

    % --- Waypoint Navigation Logic ---
    % Check if the vehicle has reached the current waypoint
    distanceToWaypoint = norm([waypoints(1, currentWaypointIndex) - trueState(1); ...
                                waypoints(2, currentWaypointIndex) - trueState(3)]);

    if distanceToWaypoint < waypointProximityThreshold(currentWaypointIndex)
        disp(['Reached waypoint: ', num2str(currentWaypointIndex)]);
        currentWaypointIndex = currentWaypointIndex + 1;
        if currentWaypointIndex > length(waypointProximityThreshold)
            disp('All waypoints reached. Ending simulation.');
            break; % Exit loop if all waypoints are done
        end
    end

    % --- 1. Simulate True System Dynamics ---
    % The 'simulateVehicleDynamics' function represents the "real world"
    nextTrueState = simulateVehicleDynamics(trueState, samplingTime, currentControlInput);

    % --- 2. Generate Noisy Measurements ---
    % Simulate GPS measurements (position)
    trueMeasurement = C_discrete * nextTrueState;
    measurementNoise = gpsNoiseMagnitude * (rand(2,1) - 0.5) * 2; % Uniform noise [-mag, mag]
    currentMeasurement_z = trueMeasurement + measurementNoise;

    % Simulate Accelerometer measurements (noisy reading of control input)
    inputNoise = accelNoiseMagnitude * (rand(2,1) - 0.5) * 2;
    currentNoisyInput_w = currentControlInput + inputNoise;

    % --- 3. Kalman Filter Estimation ---
    [nextEstimatedState, P_estimateErrorCov] = kalmanFilter( ...
        estimatedState, P_estimateErrorCov, ...
        A_discrete, B_discrete, C_discrete, ...
        R_processNoiseCov, Q_measurementNoiseCov, ...
        currentMeasurement_z, currentNoisyInput_w);

    % --- 4. Calculate Control Input ---
    % The controller uses the *estimated state* to decide on the next action.
    % (Targeting the current waypoint)
    nextControlInput = calculateControlInput(nextEstimatedState, waypoints(:, currentWaypointIndex));

    % --- Store Data for Plotting ---
    % (If not pre-allocating, grow arrays. Less efficient but simpler for fewer steps)
    allTrueStates(:, timeStep) = trueState;         % Store current true state before update
    allMeasurements(:, timeStep) = currentMeasurement_z;
    allEstimatedStates(:, timeStep) = estimatedState; % Store current estimate before update
    allControlInputs(:, timeStep) = currentControlInput;


    % --- Update States for Next Iteration ---
    trueState = nextTrueState;
    estimatedState = nextEstimatedState;
    currentControlInput = nextControlInput;


    % --- Real-time Plotting (Figure 1) ---
    if mod(timeStep, 5) == 0 || timeStep == 1 % Plot every 5 steps or first step
        plot(allTrueStates(1,1:timeStep), allTrueStates(3,1:timeStep), 'k-', 'LineWidth', 2); % Actual trajectory
        plot(currentMeasurement_z(1), currentMeasurement_z(2), 'r.', 'MarkerSize', 15);      % Current GPS reading
        plot(estimatedState(1), estimatedState(3), 'b.', 'MarkerSize', 15);                  % Current Kalman estimate
        legend('True Path', 'GPS Measurement', 'Kalman Estimate', 'Location','northwest', 'AutoUpdate','off');
        title(['Vehicle Simulation (Time: ', num2str(timeStep*samplingTime), 's)']);
        xlabel('X Position (units)');
        ylabel('Y Position (units)');
        drawnow;
    end

    if timeStep > 2000 % Safety break for very long simulations
        disp('Simulation reached maximum allowed steps.');
        break;
    end
end
disp('Simulation loop finished.');

% Ensure legend includes the final plotted points if loop ended early
if exist('allTrueStates','var')
    figure(1); % Bring to front
    h_true = plot(allTrueStates(1,1:timeStep), allTrueStates(3,1:timeStep), 'k-', 'LineWidth', 2);
    h_meas = plot(allMeasurements(1,timeStep), allMeasurements(2,timeStep), 'r.', 'MarkerSize', 15); % Last measurement
    h_est  = plot(allEstimatedStates(1,timeStep), allEstimatedStates(3,timeStep), 'b.', 'MarkerSize', 15); % Last estimate
    legend([h_true, h_meas, h_est],'True Path', 'Last GPS Meas.', 'Last Kalman Est.', 'Location','northwest');
    hold off;
end


%% Post-Simulation Results Plotting (Figure 2)
disp('Plotting results...');
if exist('allTrueStates','var') && timeStep > 1
    timeVector = (1:timeStep) * samplingTime;

    figure(2); clf;

    % --- Position Plots ---
    subplot(2,1,1);
    plot(timeVector, allMeasurements(1,:), 'r:', 'DisplayName', 'X Measured (GPS)'); hold on;
    plot(timeVector, allMeasurements(2,:), 'm:', 'DisplayName', 'Y Measured (GPS)');
    plot(timeVector, allTrueStates(1,:), 'k-', 'LineWidth', 1.5, 'DisplayName', 'X True');
    plot(timeVector, allTrueStates(3,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Y True');
    plot(timeVector, allEstimatedStates(1,:), 'b--', 'LineWidth', 1.5, 'DisplayName', 'X Estimated');
    plot(timeVector, allEstimatedStates(3,:), 'c--', 'LineWidth', 1.5, 'DisplayName', 'Y Estimated');
    title('System Positions: True, Measured, and Estimated');
    xlabel('Time (seconds)');
    ylabel('Position (units)');
    legend show;
    grid on;

    % --- Velocity Plots ---
    subplot(2,1,2);
    plot(timeVector, allTrueStates(2,:), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Vx True'); hold on;
    plot(timeVector, allTrueStates(4,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Vy True');
    plot(timeVector, allEstimatedStates(2,:), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Vx Estimated');
    plot(timeVector, allEstimatedStates(4,:), 'c--', 'LineWidth', 1.5, 'DisplayName', 'Vy Estimated');
    title('System Velocities: True and Estimated');
    xlabel('Time (seconds)');
    ylabel('Velocity (units/s)');
    legend show;
    grid on;
else
    disp('Not enough data to plot results or simulation did not run.');
end

disp('Script finished.');

%% Local Functions (formerly separate .m files)

%==========================================================================
% FUNCTION: simulateVehicleDynamics
% Purpose: Simulates the "true" physical movement of the vehicle.
% Inputs:
%   currentState_x      - Current true state [px, vx, py, vy]'
%   timeStep_T          - Sampling time
%   controlInput_u      - Applied control input (accelerations [ax; ay]')
% Output:
%   nextState_dx        - Next true state
%==========================================================================
function nextState_dx = simulateVehicleDynamics(currentState_x, timeStep_T, controlInput_u)
    % This model uses a discrete-time update based on constant acceleration
    % over the interval timeStep_T. This is effectively what A_discrete and
    % B_discrete from the main script would do if the model was perfect.
    % For teaching, this emphasizes that the filter's model (A_d, B_d)
    % is an *approximation* of reality. If this function used a more complex
    % or different model, the Kalman filter would have to cope with that mismatch.

    % State transition matrix for the true system (can be same as filter's A_d)
    A_true_model = [1 timeStep_T 0 0;
                    0 1          0 0;
                    0 0          1 timeStep_T;
                    0 0          0 1];
    % Input matrix for the true system (can be same as filter's B_d part)
    B_true_model = [0.5*timeStep_T^2  0; % More accurate: position = v*T + 0.5*a*T^2
                    timeStep_T        0;
                    0                 0.5*timeStep_T^2;
                    0                 timeStep_T];
    % If using the simplified B_d = T*B_continuous from main script:
    % B_true_model = [0          0;
    %                 timeStep_T 0;
    %                 0          0;
    %                 0          timeStep_T];

    nextState_dx = A_true_model * currentState_x + B_true_model * controlInput_u;
end

%==========================================================================
% FUNCTION: kalmanFilter
% Purpose: Implements the discrete Kalman filter algorithm.
% Inputs:
%   x_est_prev          - Previous state estimate
%   P_cov_prev          - Previous estimate error covariance
%   A_d, B_d, C_d       - Discrete system matrices
%   R_processNoise      - Process noise covariance (model uncertainty)
%   Q_measurementNoise  - Measurement noise covariance (sensor uncertainty)
%   z_measurement       - Current noisy measurement from sensors
%   w_noisyInput        - Current noisy measurement of control input
% Outputs:
%   x_est_curr          - Current (a posteriori) state estimate
%   P_cov_curr          - Current (a posteriori) estimate error covariance
%==========================================================================
function [x_est_curr, P_cov_curr] = kalmanFilter(x_est_prev, P_cov_prev, ...
                                               A_d, B_d, C_d, ...
                                               R_processNoise, Q_measurementNoise, ...
                                               z_measurement, w_noisyInput)
    % --- Prediction Step ---
    % Predict current state based on previous estimate and control input
    x_predicted = A_d * x_est_prev + B_d * w_noisyInput;
    % Predict current estimate error covariance
    P_predicted = A_d * P_cov_prev * A_d' + R_processNoise;

    % --- Update (Correction) Step ---
    % Calculate measurement residual (innovation)
    y_innovation = z_measurement - C_d * x_predicted;
    % Calculate innovation covariance
    S_innovationCov = C_d * P_predicted * C_d' + Q_measurementNoise;
    % Calculate Kalman Gain
    K_kalmanGain = P_predicted * C_d' * inv(S_innovationCov); %#ok<MINV> Standard, but pinv for robustness if S is ill-conditioned

    % Update state estimate with measurement
    x_est_curr = x_predicted + K_kalmanGain * y_innovation;
    % Update estimate error covariance
    P_cov_curr = (eye(size(A_d,1)) - K_kalmanGain * C_d) * P_predicted;
end

%==========================================================================
% FUNCTION: calculateControlInput
% Purpose: Calculates the control input (accelerations) to steer the vehicle.
% Inputs:
%   currentEstimatedState_x - Current estimated state [px, vx, py, vy]'
%   targetWaypoint_r        - Target waypoint coordinates [x_target; y_target]
% Output:
%   controlSignal_u         - Calculated control input [ax; ay]'
%==========================================================================
function controlSignal_u = calculateControlInput(currentEstimatedState_x, targetWaypoint_r)
    % Simple PD-like controller
    % Proportional term based on position error
    % Derivative-like term based on current velocity (acts as damping)

    kp = 0.05; % Proportional gain for position error
    kd = -0.5; % Gain for velocity (damping term; negative because it opposes current velocity to reduce overshoot or acts as drag)

    % Error in x and y positions
    error_px = targetWaypoint_r(1) - currentEstimatedState_x(1);
    error_py = targetWaypoint_r(2) - currentEstimatedState_x(3);

    % Current velocities
    vel_x = currentEstimatedState_x(2);
    vel_y = currentEstimatedState_x(4);

    % Calculate control signals (accelerations)
    u_x = error_px * kp + vel_x * kd;
    u_y = error_py * kp + vel_y * kd;

    controlSignal_u = [u_x; u_y];

    % Optional: Limit control input magnitude
    max_accel = 10.0; % Max acceleration in any direction (units/s^2)
    if norm(controlSignal_u) > max_accel
        controlSignal_u = controlSignal_u / norm(controlSignal_u) * max_accel;
    end
end