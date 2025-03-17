%% EGH445 - Modern Control
%--------------------------------------------------------------------------
% ModernControl
% **** Main Functions
%     **** Discretisation
%
% Author: Aaron McFadyen
%--------------------------------------------------------------------------


%% ----------------- Discretisation Examples (EGH445) --------------------
% This example demonstrates the approximation error when converting from
% continuous time to discrete time state-space systems. The exact matrix
% exponential is compared to an nth order approximation. A comparison is 
% also made with 3 different sampling times. The system is given by a 
% continuous time state-space model.

clear all; close all; clc;

% Setup Simulation Time and Sampling Time
Tm = 10;             	% Simulation Seconds
m  = 2;                 % Matrix Exp Approximation Order
n  = 10;             	% Sampling Factor (comparison purposes)
T1 = 1.0;             	% Sampling Time 1
T2 = T1/n;             	% Sampling Time 2
T3 = T1/(n/2);         	% Sampling Time 3

% Set Continuous Plant and Input Matrices
[A,B,C,D] = GetPlantModel('mass-spring');

% Find Continuous State-Space Model, Transfer Function and Eigenvalues
[num,den]   = ss2tf(A,B,C,D);
Gs          = tf(num,den)
E           = eig(Gs);

%-------------------- Find Discrete Plant and Input Matrices --------------
% Approximate Plant and Input Matrices (Mat Exp nth Order)
G1 = ApproxMatExp(A, T1, m, 1)                 % Plant (D/z) - Approximate Calculation (3rd Order)
H1 = inv(A)*(G1 - eye(size(G1)))*B      	% Plant (D/z) - Approximate Calculation (3rd Order)
E1 = eig(G1);

% Exact Plant and Input Matrices (Laplace)
syms s; 
Ga = (1-T1)*exp(-T1) + 2*T1*exp(-T1);
Gb = -T1*exp(-T1);
Gc = T1*exp(-T1);
Gd = (1-T1)*exp(-T1);
G1L = [Ga Gc;Gb Gd]
GL = inv([s 0;0 s]*eye(2)-A)

Ha = 1 - (exp(-T1)*(1+T1));
Hb = T1*exp(-T1);
H1L = [Ha;Hb]
H1LL = inv(A)*(G1L - eye(size(G1L)))*B  
Hl = GL*B./s
El = eig(GL);

% Exact Plant and Input Matrices (Matrix Transform)
[P,L] = eig(A);
if det(sym(P))~=0
    Q = exp(L.*T1).*eye(size(A));
    G1T = P*Q*inv(P);
    H1T = inv(A)*(G1T-eye(2))*B;
else
    disp('A Matrix is singular or repeated eigenvalues - invalid method')
end


% Exact Plant and Input Matrices
G2 = expm(A.*T1)                            % Input (D/z) - Exact Calculation ()
H2 = inv(A)*(G2 - eye(size(G2)))*B          % Input (D/z) - Exact Calculation ()
E2 = eig(G2);

% Sampling Time 1 (exact G and H)
G3 = expm(A.*T2)                            % Input (D/z) - Exact Calculation ()
H3 = inv(A)*(G3 - eye(size(G3)))*B          % Input (D/z) - Exact Calculation ()
E3 = eig(G3);

% Sampling Time 2 (exact G and H)
G4 = expm(A.*T3)                            % Input (D/z) - Exact Calculation () 
H4 = inv(A)*(G4 - eye(size(G4)))*B          % Input (D/z) - Exact Calculation ()
E4 = eig(G4);


C = C;
D = D;

% Find and Display System Eigenvalues
disp(horzcat('Eigenvalues (exact) Continuous System: ',num2str(E(1)),',',num2str(E(2))));
disp(horzcat('Eigenvalues (approx) Discrete System - Sampling Time ',num2str(T1),' sec:',num2str(E1(1)),',',num2str(E1(2))));
disp(horzcat('Eigenvalues (exact) Discrete System - Sampling Time ',num2str(T1),' sec:',num2str(E2(1)),',',num2str(E2(2))));
disp(horzcat('Eigenvalues (exact) Discrete System - Sampling Time ',num2str(T2),' sec:',num2str(E3(1)),',',num2str(E3(2))));
disp(horzcat('Eigenvalues (exact) Discrete System - Sampling Time ',num2str(T3),' sec:',num2str(E4(1)),',',num2str(E4(2))));

%---------- Find Discrete State-Space Model and Transfer Function ---------
% Approximate Plant and Input Matrices (nominal sampling time)
Z1          = ss(G1,H1,C,D,T1)
[num,den] = ss2tf(G1,H1,C,D,1);             % Last input defines the input number (for TF)
Gz1         = tf(num,den,T1)                % Last input defines the sampling time(for SS)

% Exact Plant and Input Matrices (Laplace)
Z1L         = ss(G1L,H1L,C,D,T1)
[num,den]   = ss2tf(G1L,H1L,C,D,1);     	% Last input defines the input number (for TF)
Gz1L        = tf(num,den,T1)                % Last input defines the sampling time(for SS)

% Exact Plant and Input Matrices (Matrix Transform)
if det(sym(P))~=0
    Z1T        = ss(G1T,H1T,C,D,T1)
    [num,den]  = ss2tf(G1T,H1T,C,D,1);     	% Last input defines the input number (for TF)
    Gz1T       = tf(num,den,T1)          	% Last input defines the sampling time(for SS)
end
% Exact Plant and Input Matrices (nominal sampling time)
Z2          = ss(G2,H2,C,D,T1)           
[num,den] = ss2tf(G2,H2,C,D,1);        
Gz2         = tf(num,den,T1)

% Sampling Time 1 (exact G and H)
Z3          = ss(G3,H3,C,D,T2)
[num,den]   = ss2tf(G3,H3,C,D,1);          
Gz3         = tf(num,den,T2)            

% Sampling Time 2 (exact G and H)
Z4          = ss(G4,H4,C,D,T3)
[num,den]= ss2tf(G4,H4,C,D,1);          
Gz4         = tf(num,den,T3)          

%------------------------- Simulate Systems -------------------------------
% Setup Time Vectors
t1 = 0:T1:Tm;
t2 = 0:T2:Tm;
t3 = 0:T3:Tm;


figure(10);hold on
impulseplot(Gs,t1(end),'k-')
%[yi,xi,ti] = impulse(Gz1,t1);
%stairs(xi,yi,'r.-')
[yi,xi,ti] = impulse(Gz1L,t1);
stairs(xi,yi,'r.-')
l=legend('$G(s)$','$G(z)$');
%l=legend('$G(s)$','$G1(z)$','$G2(z)$','$G2(z)$');
set(l,'interpreter','latex');

figure(1);hold on

% Impulse Response (unit pulse)
impulseplot(Gs,t1(end),'k-')
%impulseplot(Gz1,t,'ro-');
[yi,xi,ti] = impulse(Gz1,t1);
stairs(xi,yi,'r.-')
[yi,xi,ti] = impulse(Gz2,t1);
stairs(xi,yi,'bs-')
[yi,xi,ti] = impulse(Gz3,t2);
stairs(xi,yi,'g.-')
[yi,xi,ti] = impulse(Gz4,t3);
stairs(xi,yi,'b--')

l=legend('$G(s)$','$G1(z)$','$G2(z)$','$G3(z)$','$G4(z)$');
%l=legend('$G(s)$','$G1(z)$','$G2(z)$','$G2(z)$');
set(l,'interpreter','latex');

% Step Response (unit pulse train)
figure(2);hold on
stepplot(Gs,t1(end),'k-')
%stepplot(Gz1,t,'ro-')
[ys,xs,ts] = step(Gz1,t1);
stairs(xs,ys,'r.-')
[ys,xs,ts] = step(Gz2,t1);
stairs(xs,ys,'bs-')
[ys,xs,ts] = step(Gz3,t2);
stairs(xs,ys,'g.-')
[ys,xs,ts] = step(Gz4,t3);
stairs(xs,ys,'b--')

l=legend('$G(s)$','$G1(z)$','$G2(z)$','$G3(z)$','$G4(z)$');
%l=legend('$G(s)$','$G1(z)$','$G2(z)$','$G2(z)$');
set(l,'interpreter','latex');

% ------------------------ Simulation Resources ---------------------------
%ui = [1 zeros(1,N)];           % Input Vector (unit impulse)
%us = [1 ones(1,N)];        	% Input Vector (unit step)
%yi = filter(num2,den2,ui);
%stem(t,yi,':ob')
%ys = filter(num2,den2,us);
%stem(t,ys,':ob')

%% --------------- Difference Equation Example (EGH445) -------------------
% This example demonstrates how to simulate a difference equation
% (recursion). A simple mass-spring-damper model is used and the difference
% equations have been derived by assuking a first order (Euler Forward)
% approximation to the discrete derivative.

% Clear Variables
clear x u t

% Timing 
T  = T1/5;
Tf = 10;
t  = 0:T:Tf-T;
 
% System Parameters (mass = 1, spring = 1, damper = 2) 
a = -1; 
b = -2;
c = 1;
 
% Input
u = ones(length(t),1);           % Step Input (train)
%u = [1 zeros(length(t)-1,1)'];  % Impulse Input (train)

% Initial Conditions 
x(1,1) = 0; 
x(2,1) = 0;

% Simulate (recursion) 
for k = 2:1:length(t)
    x(1,k) = x(1,k-1) + T*x(2,k-1); 
    x(2,k) = T*a*x(1,k-1) + (1+T*b)*x(2,k-1) + T*c*u(k-1);
end

% Plot Response
figure(3); hold on
%plot(t,x(1,:),'ro--','linewidth',0.5)
%plot(t,x(2,:),'bo--','linewidth',0.5)
stem(t,x(1,:),'ro-','linewidth',0.5)
stem(t,x(2,:),'bo-','linewidth',0.5)
l=legend('$x_1(t)$','$x_2(t)$');
set(l,'interpreter','latex');
xlabel('Time (s)'); ylabel('States (m)');


figure(4); hold on; grid on
plot(x(1,:),x(2,:),'r.-','linewidth',2)
plot(x(1,end),x(2,end),'ro','linewidth',2)
plot(x(1,1),x(2,1),'r*','linewidth',2)
xlabel('State 1 (m)'); ylabel('State 2 (m)');
