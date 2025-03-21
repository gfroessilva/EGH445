%% EGH445 - Modern Control
%--------------------------------------------------------------------------
% Digital Control
% **** Main Functions
%     **** System Response, Stability (and Observability, Controllability)
%
% Author: Aaron McFadyen
%--------------------------------------------------------------------------


%% System Analysis - Pulse Transfer Funtion, Stability (requires update)
clear all;
close all;
clc;

% Useful Symbols
syms s z 

% Sampling Time
T = [1.0 0.5 0.1];  % Sampling Period (s) (can define many)
f = 1./T;           % Sampling Frequency (Hz)
Tm = 10;
t = 0:T:Tm;

% Useful Symbolic Variables
syms s z

% System Model Parameters
A = [0 1;-1 -2];    B = [0;1];          % Model Parameters (Input Eqn)
C = [0 1];          D = 0;              % Model Parameters (Output Eqn)
[A,B,C,D] = GetPlantModel('mass-spring');   

% Feedback Model (unit/direct output feedback)
Hs.tf = tf([0 1],[0 1])
Hz.tf = tf([0 1],[0 1])


%% Analyse System(s) - Poles (s-z plane/mapping)
CSS = ss(A,B,C,D);                      % Continuous State Space Models
[num,den] = ss2tf(A,B,C,D);
Gs.tf = tf(num,den);                    % Continuous Transfer Function
Gs.poles = roots(den);                  % Continuous System Poles
for k=1:1:length(T)
    % Recall z = e^(Ts)
    temppoles(k,:) = exp(T(k).*ones(size(Gs.poles)).*Gs.poles)'; % s -> z mapping
end
disp(num2str(temppoles))

%% Analyse System(s) - Stability, Response (open and closed)
%** Also finds Controllability and Observability (see commented out) **

for k=1:1:1%length(T) % Loop through/evaluate different sampling times
    % Timing
    t = 0:T(k):Tm;
    
    % -------------- System Modelling - State Space -----------------------
    % Get Continuous Model
    Gs.ss = ss(A,B,C,D)                       % Continuous State Space Models
    [num,den] = ss2tf(A,B,C,D);
    Gs.tf = tf(num,den);                    % Continuous Transfer Function
    Gs.poles = roots(den) 
    figure(k*10); hold on; pzmap(Gs.tf);
    grid on; title('s Plane','interpreter','latex')
    %ExportFigJpg(horzcat('splane_T',num2str(k)));
    
    % Get Continuous Model -> Discrete (option 1)
    G   = expm(A.*T(k))                     % Discrete State Mtx - Exact Calculation                      
    H   = inv(A)*(G - eye(size(G)))*B       % Discrete Input Mtx - Exact Calculation 
    Fz.ss = ss(G,H,C,D,T(k));               % Discrete State Space Models   
    
    % Get Continuous Model -> Discrete (option 2)             
    [G,H] = c2d(A,B,T(k));                  % Discrete State and Input Mtx
    Fz.ss   = ss(G,H,C,D,T(k));             % Discrete State Space Model
    
    %--------------------
    % Find (open Loop) Pulse Transfer Function (option 1)
    [num,den] = ss2tf(G,H,C,D);
    Fz.tf     = tf(num,den,T(k))
    Fz.poles  = roots(den)
    Fz.zeros  = roots(num)
    %Fz.poles  = roots(poly(G));
    
    % Find (Closed Loop) Pulse Transfer Function (option 1)
    Hz.tf           = tf([0 1],[0 1],T(k))
    Fz.tfclosed     = feedback(Fz.tf,Hz.tf)
    Fz.polesclosed  = roots(cell2mat(Fz.tfclosed.den))
    
    figure(k*10 + 1); hold on; 
    h = zplane(Fz.zeros,Fz.poles);       
    grid on; title('z Plane',...
        'interpreter','latex')
    %ExportFigJpg(horzcat('zplane_T',num2str(k)));
    
    % Find Pulse Transfer Function (option 2)
    F1 = simplify(vpa(C*inv([z 0;0 z]*eye(length(G))-G)*H + D,2))
       
    % Approximate Transfer Function (option 3 - check manual calculations)
    Ghat = round(G,2);                          % Round G to 2 d.p
    Hhat = round(H,2);                          % Round H to 2 d.p
    F2 = simplify(vpa(C*inv([z 0;0 z]*eye(length(Ghat))-Ghat)*Hhat + D,2))
    Gdash = [z 0;0 z]*eye(length(Ghat))-Ghat;   % (zI-G)
    Gadj  = vpa(adjoint(Gdash),2);              % Adjoint    (zI-G)
    Gdet  = vpa(det(Gdash),2);                  % Determinant(zI-G)
    Ginv  = vpa(inv(Gdash),2);                  % Inverse    (zI-G)
    
    Fnum = vpa(C*Gadj*Hhat,2);                  % Numerator     F(z)
    Fnum = double(coeffs(Fnum,z));    	
    Fden = vpa(Gdet,2);                         % Denominator   F(z)
    Fden = double(coeffs(Fden,z));    
    F3 = tf(fliplr(Fnum),fliplr(Fden),T(k))
    
    Gk = vpa(Ginv*z,2)
    Gint = vpa(Gk*Hhat,2)
    
	% ------------------------- Find Stability ----------------------------
    % Stability (Open Loop)
    E = eig(G);                 % Check Poles (Eigenvalues G)
    disp(horzcat('Open Loop Poles: ', num2str(E'), ' (K = ', num2str(0),')'))
    
    % Stability (Closed Loop)
    syms k1 k2
    K.sym = [k1 k2];
    Ec.sym = vpa(eig(G-H*K.sym),2);

    k1 = 0;k2=1; K.val = [k1 k2]; 
    Ec.val = eval(subs(Ec.sym));
    disp(horzcat('Closed Loop Poles: ', num2str(Ec.val'),' (K = ', num2str(K.val),')')); 

    figure(99); hold on
    stepplot(Gs.tf,t(end),'k-')
    %stepplot(Gz1,t,'ro-')
    [ys,xs,ts] = step(Fz.tf,t);
    stairs(xs,ys,'r-')
    [ys,xs,ts] = step(F3,t);
    stairs(xs,ys,'b--')    

    % ----------------- Controllability and Observability -----------------
%     % Controllability 
%     Cgh.mtx = [H G*H];
%     Cgh.det = det(sym(Cgh.mtx));
%     Cgh.rnk = rank(Cgh.mtx);
%     %Cgh.rnk = rank(ctrb(G,H));
%     if Cgh.det == 0
%         disp(horzcat('System Uncontrollable (Rank: ',num2str(Cgh.rnk),')'))
%     else
%         disp(horzcat('System Controllable (Rank: ',num2str(Cgh.rnk),')'))
%     end
%     
%     % Observability
%     Ogh.mtx = vertcat(C, C*G); 	% This is the same as lecture slides
%     Ogh.mtx = [C; (C*G)];
%     Ogh.mtx = [C' (C*G)']';    	% This is the alternate formulation
%     Ogh.det = det(sym(Ogh.mtx));
%     Ogh.rnk = rank(Ogh.mtx);
%     %Ogh.rnk = rank(obsv(G,C));
%     if Ogh.det == 0
%         disp(horzcat('System Unobservable (Rank: ',num2str(Ogh.rnk),')'))
%     else
%         disp(horzcat('System Observable (Rank: ',num2str(Ogh.rnk),')'))
%     end    
%     disp('')
     
    
end




%% System Response
%clear all
%close all
clc

% Timing 
T  = 1;
Tf = 10;
t  = 0:T:Tf;
 
% System Parameters (mass = 1, spring = 1, damper = 2) 
c = 0.3679; 
a = -log(c)/T;

% Input
u = ones(length(t),1);       % Step Input (train)

% Initial Conditions 
y(1,1) = 0; 

% Simulate (recursion) 
for k = 2:1:length(t)
    y(1,k) = ((k-1)*T)*exp(-a*(k-1)*T); 
end

figure(999);hold on
stairs(t,y,'-','Color',[0.8 0.8 0.8],'LineWidth',3)
stepplot(Fz.tf,t(end),'b--')
stepplot(Gs.tf,t(end),'k-')
l=legend('$y_1(kT)$','$y_2(kT)$','$y(t)$');
set(l,'interpreter','latex');

return


%% Stability Resources
% zplane(b,a) 
% [vz,vp,vk] = zplane(d) 
% sys = tf([0.04798 0.0464],[1 -1.81 0.9048],0.1);
% P = pole(sys) 

% http://control.dii.unisi.it/sdc/altro/TabellaTrasformataZ.pdf - Tables

%% Control and Observeability Resources (Analysis and Controller Design)
% Find controllability matrix of the state-space LTI system
%Co     = ctrb(A,B)                 % Option 1
%Co     = ctrb(sys)                 % Option 2
%Co     = ctrb(sys.A,sys.B);        % Find controllability matrix of cascaded system (series)
%unco    = length(A) - rank(Co)     % Determine the number of uncontrollable states.

% Find observability matrix of the state-space LTI system
%Ob      = obsv(A,C)                % Option 1
%Ob      = obsv(sys)                % Option 2
%unob    = length(A)-rank(Ob)       % Determine the number of unobservable states.


%% System Response Resources
%[h,t] = impz(b,a)
%[h,t] = stepz(b,a)
