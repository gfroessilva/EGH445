%% EGH445 - Modern Control
%--------------------------------------------------------------------------
% Digital Control
% **** Main Functions
%       **** Direct/Discrete Design and Emulation Comparisons
%     
% Author:   Dr Aaron McFadyen
%--------------------------------------------------------------------------



clear all; close all; clc;


%% Parameters
% Setup Simulation Time and Sampling Time
T = [1.0 0.5 0.1];  % Sampling Period    (s)
%T = 1;
f = 1./T;           % Sampling Frequency (Hz)
Tm = 20;            % Simulation Time    (s)
t = 0:T:Tm;         % Simulation Samples (s)

% Useful Symbolic Variables
syms s z k0 k1

% Useful Plotting Tool
colors.full  = {'b','r','g'};
colors.light = [204 229 255;255 204 204;204 255 204]./255; 

%% System Modelling, Analysis and Control Design
% For different sampling times, model the system (continuous and discrete)
% then analyse the system behaviour (stability and reponse). Deisgn
% controllers for the discrete time system using direct and emulation
% design approaches. Compare results.
                                        
for k=1:1:1%length(T)
    
    % Timing and Parameters
    t = 0:T(k):Tm;
    
    
    % ************************ System Modelling **************************
    % System Model Parameters
    % A = [0 1;-1 -2];    B = [0;1];            % Model Parameters (Input Eqn)
    % C = [0 1];          D = 0;                % Model Parameters (Output Eqn)
    
    [A,B,C,D] = GetPlantModel('mass-spring');   % Select Model  
    C = [1 0];                                  % Change Output Equation 
                                                % (measure distance not velocity)

    % System Modelling - Continuous Time 
    % Get Continuous System State Space Model
    Gs.ss = ss(A,B,C,D)                         % Continuous State Space Models
    
    % Get Continuous (Open Loop) Pulse Transfer Function
    [num,den] = ss2tf(A,B,C,D);
    Gs.tf = tf(num,den);                        % Continuous Transfer Function
   
    
    % System Modelling - Discrete Time
    % Get Discrete System State Space Model (option 1)
    G       = expm(A.*T(k))                   	% Discrete State Mtx - Exact Calculation                      
    H       = inv(A)*(G - eye(size(G)))*B       % Discrete Input Mtx - Exact Calculation 
    Fz.ss   = ss(G,H,C,D,T(k));                 % Discrete State Space Models   
    
    % Get Discrete System State Space Model (option 2)             
    [G,H]   = c2d(A,B,T(k));                    % Discrete State and Input Mtx
    Fz.ss   = ss(G,H,C,D,T(k));                 % Discrete State Space Model
    
    % Get Discrete (Open Loop) Pulse Transfer Function
    [num,den] = ss2tf(G,H,C,D);
    Fz.tf     = tf(num,den,T(k))
    %F1       = simplify(vpa(C*inv([z 0;0 z]*eye(length(G))-G)*H + D,2))% Pulse Transfer Function (option 2)
       
    % ************************ System Analysis **************************
    % System Poles - Continuous Time 
    den = cell2mat(Gs.tf.den);
    num = cell2mat(Gs.tf.num);
    Gs.poles = roots(den);
    %Gs.poles = roots(poly(A));
    Gs.zeros = roots(num); 
    
    % Stability Analysis (Open Loop)
    Gs.eig = eig(A);                 % Check Poles (Eigenvalues G)
    for i=1:1:length(Gs.eig)
        disp(horzcat('Continuous Open Loop Pole (Eigenvalue) ',num2str(i),': ', num2str(Gs.poles(i)),'(',num2str(Gs.eig(i)),')'))
    end
    disp(' ')
    
	figure(k*10); hold on; grid on; 
    pzmap(Gs.tf);
    title('s Plane','interpreter','latex')
    %ExportFigJpg(horzcat('splane_T',num2str(k)));
    
    % System Poles - Discrete Time
    % Find (Closed Loop) Pulse Transfer Function (option 1)
    %[mapping] = GetPoleMapping(A,B,C,D,T(k));
    %Fz.poles = mapping.poles;
    den = cell2mat(Fz.tf.den);
    num = cell2mat(Fz.tf.num);
    Fz.poles  = roots(den);
    %Fz.poles = roots(poly(G));
    Fz.zeros  = roots(num);
        
    % Stability Analysis (Open Loop)
    Fz.eig = eig(G);                 % Check Poles (Eigenvalues G)
    for i=1:1:length(Fz.eig)
        disp(horzcat('Discrete   Open Loop Pole (Eigenvalue) ',num2str(i),': ', num2str(Fz.poles(i)),'(',num2str(Fz.eig(i)),')'))
    end
    disp(' ')
    
	figure(k*10 + 1); hold on; grid on;
    h = zplane(Fz.zeros,Fz.poles);       
    title('z Plane', 'interpreter','latex')
    %ExportFigJpg(horzcat('zplane_T',num2str(k)));
    
    % Controllability 
    Cgh.mtx = ctrb(G,H);
    Cgh.det = det(sym(Cgh.mtx));
    Cgh.rnk = rank(Cgh.mtx);
    if Cgh.det == 0
        disp(horzcat('System Uncontrollable (Rank: ',num2str(Cgh.rnk),')'))
        disp('Control Design Incomplete - Exiting Design Procedure')
        return
    else
        disp(horzcat('System Controllable (Rank: ',num2str(Cgh.rnk),')'))
    end
    
    %  Observability
    Ogh.mtx = obsv(G,C);
    Ogh.det = det(sym(Ogh.mtx));
    Ogh.rnk = rank(Ogh.mtx);
    if Ogh.det == 0
        disp(horzcat('System Unobservable (Rank: ',num2str(Ogh.rnk),')'))
        disp('Control Design Incomplete - Exiting Design Procedure')
        return
    else
        disp(horzcat('System Observable   (Rank: ',num2str(Ogh.rnk),')'))
    end
    disp(' ')


    % Get Feedback Transfer Functions 
    % Consider unity gain feedback (for reference)
    Hz.tf = tf([0 1],[0 1],T(k));   % Discrete
    Hs.tf = tf([0 1],[0 1]);        % Continuous  
    
    
    % **************** System Control Design - Emulation ******************
    Kc = 2.8; % Continuous Controller Gain
    %Kc = 1;
    
    % Define Continuous Controller Types
    %Gc.gain = tf([0 Kc],1);
    Gc.lead = tf(Kc.*[1 2],[1 4]);
    %Gc.lag = tf(Kc/2.*[1 10],[1 1]);
    %Gc.leadlag = series(Gc.lead,Gc.lag);
    %Gc.pd = tf(pid(0.5,0.0,2.0,0));    
    
    % Find Discrete Controller
    Gc.tf = Gc.lead; 
    methods = {'tustin','zoh','matched'};% Define Emulation Methods
    [F] = Emulation(Gc,Gs,Fz,t,methods,colors);
    
    
    % ***************** System Control Design - Direct ********************
    % Performance Metrics
    wn   = 2;
    zeta = 0.8;
    
    % Find Desired Char Equation
    CharEqn.desired.s   = [1 2*zeta*wn wn^2];
    Roots.desired.s     = roots(CharEqn.desired.s);
    Roots.desired.z     = exp(roots(CharEqn.desired.s).*T(k));
    CharEqn.desired.z   = poly(Roots.desired.z);
    
    % Transform System 
    CharEqn.open.z  = poly(G);                      % Find Actual Char Equation (Open Loop)
    Gtilde  = [0 1; -flip(CharEqn.open.z(2:end))];  % G - Phase Variable/Companion
    Htilde  = [0 1]';                               % H - Phase Variable/Companion
    
    % Get Transformtion Matrix
    Cgh.orig  = ctrb(G,H);
    Cgh.tilde = ctrb(Gtilde,Htilde);
    P = Cgh.orig*inv(Cgh.tilde);              	
    
    Gp = inv(P)*G*P;
    Hp = inv(P)*H;
    Cp = C*P;
    Dp = D;
    Fz.ssP = ss(Gp,Hp,Cp,Dp,T(k));
    
    % Find Gain (Option 1 - Transforms)
    Kp = Gp(end,:)  +  flip(CharEqn.desired.z(2:end));   % Transformed System Gain
    K  = Kp*inv(P);                                      % Original System Gain    
    
    % Find Gain (Option 2 - Equating Coeffs + Sim. Eqns.)
    syms k0 k1
    Ke = [k0 k1];    
    Gk = vpa(G-H*Ke,2);
    Gkdet = vpa(det(z.*eye(length(G))-Gk),2);
    
    % Find Gain (Option 3 - Ackerman)
    Ka = place(G,H,Roots.desired.z);

    % Feed Forward Gain Adjustment

    % Continuous Feedback Control (unity feedback, discrete)
    Gs.tfclosed     = feedback(series(Gc.tf,Gs.tf),Hs.tf);
    Gs.polesclosed  = roots(cell2mat(Gs.tfclosed.den));
    Gs.stepinfo     = stepinfo(Gs.tfclosed);

    [num, den]= ss2tf(G-(H*K),H,C,D);
    Fz.tfsf = tf(num,den,T(k));
    Ks = evalfr(Gs.tfclosed,0);
    Kr = Ks/evalfr(Fz.tfsf,1);
    if T(k)==1              % Sampling Time T = 1
        KK = [1 1.05 1.1 1.15 1.2 1.25 1.3];
    elseif T(k) == 0.1      % Sampling Time T = 0.1
        KK = [1 1.6 1.625 1.65 2];    
    else
        KK = 1;
    end

    Fz.ss = ss(G-(H*K),H,C,D,T(k));
    Fz.ssff = ss(G-(H*K),Kr*H,C,D,T(k));
    
    figure(1001); hold on
    for j=1:1:length(KK)
        [x2,t2] = step(series(KK(j),Fz.ss),t);
        stairs(t2,x2,'b-','LineWidth',1);
    end
    [x1,t1] = step(Fz.ssff,t);              % Feedforward Gain = Kr
    %[x1,t1] = step(series(Kr,Fz.ss),t); 	% Feedforward Gain = Kr 
	stairs(t1,x1,'b-','LineWidth',3);
    [x2,t2] = step(series(1,Fz.ss),t);  	% Feedforward Gain = 1
	stairs(t2,x2,'LineStyle','-','Color',colors.light(1,:),'LineWidth',3);
    step(Gs.tfclosed,t,'k-');
    xlabel('Time (s)'); ylabel('Amplitude')
    
    %figure(505)
    %stairs(t2,x2,'LineStyle','-','Color',colors.light(1,:),'LineWidth',3);
    %stairs(t1,x1,'b-','LineWidth',3);
    
    % Check Closed Loop Poles Match     
    disp(horzcat('Desired System CharEqn:     ',num2str(CharEqn.desired.z)))
    disp(horzcat('Original System CharEqn:    ',num2str(poly(G-(H*K)))))
    disp(horzcat('Original System CharEqn:    ',num2str(poly(G-(H*Ka)))))
    disp(horzcat('Transformed System CharEqn: ',num2str(poly(Gp-(Hp*Kp)))))
    as=1;
            
end


%% Control Design Resources 1
% Desired Pole Location

% poles.desired.s = [-4 + 15i -4 - 15i];  % Desired Poles (s domain)
% poles.desired.z = exp(poles.desired.s.*T);      % Desired Poles (z domain)
% chareqn.desired = conv([1 poles.desired.z(1)],[1 poles.desired.z(2)])
% 
% 
% G = zpk([],[-5 -5 -10],100);
% C1 = pid(2.9,7.1);
% CL1 = feedback(G*C1,1);
% C2 = pid(29,7.1);
% CL2 = feedback(G*C2,1);
% pzplot(CL1,CL2)
% grid

%% Control Design Resources 2
% clc
% close all
% clear all
% 
% syms s z k1 k2 k3 k4
% 
% G = [0.94 0.065;-1.2 0.35];
% H = [0.035;0.066];
% K = [];
% 
% poles.actual.z  = eig(G);
% poles.desired.z = [0.5+0.2i ; 0.5-0.2i];
% 
% chareqn.open    = det(z.*eye(size(G)) - G);                              % Open   Loop Characterisitic Equation
% chareqn.closed  = det(z.*eye(size(G)) - (G - H*([k1 k2])));              % Closed Loop Characterisitic Equation
% chareqn.desired = conv([1 -poles.desired.z(1)],[1 -poles.desired.z(2)]); % Desired     Characterisitic Equation
% chareqn.desired = (z - poles.desired.z(1))*(z - poles.desired.z(2));
% 
% vpa(chareqn.open,4)
% vpa(chareqn.closed,4)
% vpa(chareqn.desired,4); simplify(vpa(chareqn.desired,4))

%% Control Design Resources 3
% G = tf([2 5 1],[1 2 3],'inputname','torque',...
%                          'outputname','velocity');
% H = zpk(-2,-10,5)
% Cloop = feedback(G,K)

% K = place(A,B,p)
% [K,prec,message] = place(A,B,p)
% 
% p = [-1 -1.23 -5.0]; % closed-loop poles at p = [-1 -1.23 -5.0]
% K = place(a,b,p)    % state-space system (a,b,c,d) with two inputs, three outputs, and three states. 
% 
% y = filter(b,a,x) %filters the input data x using a rational transfer function defined by the numerator and denominator coefficients b and a.
% 
% % Therefore, a(1) must be nonzero.

%% Control Resources  4
% Reduce System Order / Minimal Realization (pole-zero cancellation)
% sysr    = minreal(sys) %eliminates uncontrollable or unobservable state in state-space models, or cancels pole-zero pairs
% %in transfer functions or zero-pole-gain models. The output sysr has minimal order and the same response characteristics as the original model sys
% g = zpk([],1,1);
% h = tf([2 1],[1 0]);
% cloop = inv(1+g*h) * g
% cloopmin = minreal(cloop)









