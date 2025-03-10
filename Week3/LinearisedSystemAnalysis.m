
%% EGH445 Modern Control
% Nonlinear Systems and Linearisation - Linearised System Analysis

% Consider a single link robot arm system.
% Plot the phase portrait/plot for the unforced (no input) and forced
% (input) system.

clear all; close all; clc
addpath('PhasePlot')

%% System and Simulation Parameters
global m b l g

m=0.5;
b=0.3; 
l=0.4;
g=9.81;

% Simulation Time and Options
tspan = [0 15];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4); % Solver/Integrator Options

%% Initial Conditions
x0 = [-pi/2 1/(2*pi)];
ut = linspace(0,tspan(end),25);
%u = (g/l)*sin(pi/4).*ones(numel(ut),1);    % Constant input (for EP x = [pi/4, 0])
u = 0.*ones(numel(ut),1);                   % Zero input

%% Simulate System(s)
flag = 0;
[T,x]=ode45(@(t,x) robotarm(t,x,ut,u,flag), tspan, x0,opts);    % Linear (Robot Arm)

flag = 1;
[T1,x1]=ode45(@(t,x) robotarm(t,x,ut,u,flag), tspan, x0,opts);  % Nonlinear (Robot Arm)

x1=rad2deg(x1);
x = rad2deg(x);

% PLot System Responses
figure(1);clf;
title('Time Response'); hold on; grid on
subplot(211); hold on; grid on
title('Time Response');
plot(T,x(:,1),'b--');
plot(T1,x1(:,1),'r--');
legend('Linear','Nonlinear')
ylabel('q (deg)','Interpreter','Latex','FontSize',8)
xlabel('Time','FontSize',8)

subplot(212); hold on; grid on
title('Time Response');
plot(T,x(:,2),'b-');
plot(T1,x1(:,2),'r-');
legend('Linear','Nonlinear')
ylabel('$\dot{q}$ (deg/s)','Interpreter','Latex','FontSize',8)
xlabel('Time','FontSize',8)
%text(4,0.3,['k 2=',num2str(k2)],'FontSize',8

% Plot State Trajectories 
figure(2);clf;
title('State Trajectory'); hold on; grid on
%subplot(211); hold on
plot(x(:,1),x(:,2),'b--');
plot(x1(:,1),x1(:,2),'r--');
legend('Linear','Nonlinear')
xlabel('$q$ (deg)','Interpreter','Latex','FontSize',8)
ylabel('$\dot{q}$ (deg/s)','Interpreter','Latex','FontSize',8)

%% Analyse System(s)

A1 = [0 1;-g/l -b/(m*l^2)];                 e1 = eig(A1)
A2 = [0 1; g/l -b/(m*l^2)];                 e2 = eig(A2)
A3 = [0 1;-(g/l)*(sqrt(2)/2) -b/(m*l^2)];   e3 = eig(A3)

% Consider u = 0 (no input)
u = 0.*ones(numel(ut),1);               % Zero input
odefun = @(t,x) [robotarm(t,x,ut,u,1)];
figure(3); hold on;
plotpp(odefun,'tspan', tspan(end),...
    'quivercolor', [0.6,0.6,0.6],'linecolor', [0.3,0.3,0.3])
xlabel('$q$ (rad)','Interpreter','Latex','FontSize',8)
ylabel('$\dot{q}$ (rad/s)','Interpreter','Latex','FontSize',8)
title('Phase Portrait (Plot) - No Input/Unforced')

% Consider u != 0 (input)
u = (g/l)*sin(pi/4).*ones(numel(ut),1); % Constant input (for EP x = [pi/4, 0])
odefun = @(t,x) [robotarm(t,x,ut,u,1)];
figure(4); hold on;
plotpp(odefun,'tspan', tspan(end),...
    'quivercolor', [0.6,0.6,0.6],'linecolor', [0.3,0.3,0.3])
xlabel('$q$ (rad)','Interpreter','Latex','FontSize',8)
ylabel('$\dot{q}$ (rad/s)','Interpreter','Latex','FontSize',8)
title('Phase Portrait (Plot) - Input/Forced')

%% Define Robot Arm Dynamics 
function [xdot] = robotarm(t,x,ut,u,flag)

global m b l g

Tq = interp1(ut,u,t); % Evaluate control at time t

if flag % Nonlinear Model
    xdot(1) = x(2);
    xdot(2) = -(g/l)*sin(x(1))-(b/(m*l^2)*x(2)) + Tq ;
    xdot = xdot';
else % Linear Model (Equilbrium Point 3)
    xdot(1) = x(2);
    xdot(2) = -(g/l)*(sqrt(2)/2)*x(1)-(b/(m*l^2)*x(2)) + Tq ;
    xdot = xdot';
end

end
