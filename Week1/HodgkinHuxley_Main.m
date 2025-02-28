%% Clear Workspace and Close Figures
clear;
close all;

% Simulate Step Response of model
stoptime = 100; % Length of the sim [ms]
Open = sim('HodgkinHuxley_StepResponse.slx', 'Solver', 'ode15s', 'StopTime', 'stoptime');

% Plot
figure(1);
title('Open-loop response of Neuron under Step Disturbance')
plot(Open.V.Time, Open.V.Data);
ylabel('Membrane Voltage [mV]');
xlabel('Time [ms]')