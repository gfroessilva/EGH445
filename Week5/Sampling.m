%% EGH445 - Modern Control
%--------------------------------------------------------------------------
% Digital Control
% **** Main Function
%     **** Sampling and Reconstruction Examples
%         
%         Sample simple and complex signal with or without noise at frequencies at,
%         below and above the Nyquist rate. Attempt to reconstruct the sampled
%         ignals using linear, sinusoidal and optimised fitting techniques.
%
% Author: Aaron McFadyen
%--------------------------------------------------------------------------


%% Signal Construction

close all
clear all
clc

% Signal Parameters
f1  = 25;             	% Signal 1 Frequency 
f2  = 45;             	% Signal 2 Frequency 
a1  = 1;               	% Signal 1 Amplitude
a2  = 0.5;           	% Signal 2 Amplitude
fw  = 500;           	% Noise Signal Frequency
w   = 0.2;             	% Noise Signal Amplitude
recon = 'spline';
noiseon = 1;

% Timing and Signals
% Timing
ft = 1000;                  % High Sampling Frequency 
Tm = 0.1;                   % Max Simulation Time
signal.t  = 0:1/ft:Tm;      % Time Vector

% Signals
signalone.f  = @(f,t,a) a*cos(2*pi*f*t);
signaltwo.f  = @(f,t,a) a*sin(2*pi*f*t);
noise.f      = @(f,t,a) a*sin(2*pi*f*t);

signalone.s  = signalone.f(f1,signal.t,a1);      	% Signal 1
signaltwo.s  = signaltwo.f(f2,signal.t,a2);         % Signal 2
signal.s     = signalone.s + signaltwo.s;           % Combined Signal (1 + 2)

if noiseon
    %noise.sgn    = randi([-1,1],length(signal.t),1)';
    %noise.s      = noise.sgn.*(w*cos(2*pi*fw*signal.t - pi/4)); 	% Noise Signal
    %noise.s      = 0.1.*randn(length(signal.t),1) + 0;           	% Noise Signal
    signal.snom  = signalone.s + signaltwo.s;
    noise.s      = w*cos(2*pi*fw*signal.t);                       	% Noise Signal
    signal.s     = signalone.s + signaltwo.s + noise.s;             % Combined Signal (1 + 2 + Noise)
end
% Signal Labels/Names
signalone.name  = horzcat('$x_1(t)$');
signaltwo.name  = horzcat('$x_2(t)$');
noise.name      = '$w(t)$';
signal.name     = horzcat('$x(t)$');
% signalone.name = horzcat('$',num2str(a1),'\cos(2\pi',num2str(f1),'t)$');
% signaltwo.name = horzcat('$',num2str(a2),'\sin(2\pi',num2str(f2),'t)$');
% if noiseon
%     signal.name = horzcat('$x(t)$', '+', noise.name);
% else
%     signal.name = horzcat('$x(t)$');
% end

% Plot Signals and Components
figure(1);hold on;
plot(signal.t, signal.s, '-','LineWidth',4,'Color',[0.8 0.8 0.8]);
plot(signal.t, signalone.s, '--','LineWidth',2,'Color',[0.8 0.8 0.8 0.6]);
plot(signal.t, signaltwo.s, ':','LineWidth',2,'Color',[0.8 0.8 0.8 0.6]);
if noiseon
    plot(signal.t, noise.s, '-','LineWidth',1,'Color',[0.8 0.8 0.8 0.6]);
    plot(signal.t, signal.snom, '-','LineWidth',1,'Color',[0 0 0 0.6]);
    l=legend(signal.name,signalone.name,signaltwo.name,noise.name,'$\bar{x}(t)$');
else
    l=legend(signal.name,signalone.name,signaltwo.name);
end
set(l,'interpreter','latex');
xlabel('Time (s)');
ylabel('Signal Value')

%% --------------------- Complex Sampling Example -------------------------
close all
% Sampled Signal
fm = max([f1 f2]);
nyquist = 2*fm;

% Under Sampled
sampled.f1 = nyquist*(1/3);
sampled.t1 = 0:1/sampled.f1:Tm;
if noiseon
    sampled.s1 = signalone.f(f1,sampled.t1,a1) + signaltwo.f(f2,sampled.t1,a2)...
        + noise.f(fw,sampled.t1,w);
else
    sampled.s1 = signalone.f(f1,sampled.t1,a1) + signaltwo.f(f2,sampled.t1,a2);
end
[sampled.pt1,sampled.ps1] = FitFunction(sampled.t1, sampled.s1);
%sampled.p1 = interp1(sampled.t1,sampled.s1,signal.t,recon,0);

% Over Sampled (very fast)
sampled.f2 = nyquist*5;
sampled.t2 = 0:1/sampled.f2:Tm;
if noiseon
    sampled.s2 = signalone.f(f1,sampled.t2,a1) + signaltwo.f(f2,sampled.t2,a2)...
        + noise.f(fw,sampled.t2,w);
else
    sampled.s2 = signalone.f(f1,sampled.t2,a1) + signaltwo.f(f2,sampled.t2,a2);
end
sampled.p2 = interp1(sampled.t2,sampled.s2,signal.t,recon,0);

% Over Sampled (fast)
sampled.f3 = nyquist*2;
sampled.t3 = 0:1/sampled.f3:Tm;
if noiseon
    sampled.s3 = signalone.f(f1,sampled.t3,a1) + signaltwo.f(f2,sampled.t3,a2)...
        + noise.f(fw,sampled.t3,w);
else
    sampled.s3 = signalone.f(f1,sampled.t3,a1) + signaltwo.f(f2,sampled.t3,a2);
end
sampled.p3 = interp1(sampled.t3,sampled.s3,signal.t,recon,0);

% Nyquist Sampled
sampled.fn = nyquist*1;
sampled.tn = 0:1/sampled.fn:Tm;
if noiseon
    sampled.sn = signalone.f(f1,sampled.tn,a1) + signaltwo.f(f2,sampled.tn,a2)...
        + noise.f(fw,sampled.tn,w);
else
    sampled.sn = signalone.f(f1,sampled.tn,a1) + signaltwo.f(f2,sampled.tn,a2);
end
[sampled.ptn,sampled.psn] = FitFunction(sampled.tn, sampled.sn);
%p = polyfit(sampled.tn, sampled.sn,2);
%sampled.pn = polyval(p,sampled.tn);
%sampled.psn = interp1(sampled.tn,sampled.sn,signal.t,recon,0);
%sampled.ptn = sampled.tn;

figure(11);hold on;
plot(signal.t, signal.s, '-','LineWidth',4,'Color',[0.8 0.8 0.8 0.4]);
stem(sampled.t1, sampled.s1, '-o','LineWidth',1,'Color',[1 0 0 0.5]);
plot(sampled.t1, sampled.s1,'-','LineWidth',3,'Color',[1 0 0 1])
plot(sampled.pt1, sampled.ps1,'--','LineWidth',1,'Color',[1 0 0 1])
if noiseon
    plot(signal.t, signal.snom, '-','LineWidth',1,'Color',[0 0 0 0.6]);
    l=legend('$x(t)$','$x(kT)$','$\hat{x}(t)$','$\tilde{x}(t)$','$\bar{x}(t)$');
else
    l=legend('$x(t)$','$x(kT)$','$\hat{x}(t)$','$\tilde{x}(t)$');
end
set(l,'interpreter','latex');
xlabel('Time (s)');
ylabel('Signal Value')

figure(111);hold on;
plot(signal.t, signal.s, '-','LineWidth',4,'Color',[0.8 0.8 0.8 0.4]);
stem(sampled.tn, sampled.sn, '-o','LineWidth',1,'Color',[0 1 0 0.5]);
plot(sampled.tn, sampled.sn,'-','LineWidth',3,'Color',[0 1 0 1])
plot(sampled.ptn, sampled.psn,'--','LineWidth',1,'Color',[0 1 0 1])
if noiseon
    plot(signal.t, signal.snom, '-','LineWidth',1,'Color',[0 0 0 0.6]);
    l=legend('$x(t)$','$x(kT)$','$\hat{x}(t)$','$\tilde{x}(t)$','$\bar{x}(t)$');
else
    l=legend('$x(t)$','$x(kT)$','$\hat{x}(t)$','$\tilde{x}(t)$');
end
set(l,'interpreter','latex');
xlabel('Time (s)');
ylabel('Signal Value')

figure(1111);hold on;
plot(signal.t, signal.s, '-','LineWidth',4,'Color',[0.8 0.8 0.8 0.4]);
stem(sampled.t2, sampled.s2, '-o','LineWidth',1,'Color',[0 0 1 0.5]);
plot(sampled.t2, sampled.s2,'-','LineWidth',3,'Color',[0 0 1 1])
plot(signal.t, sampled.p2,'--','LineWidth',1,'Color',[0 0 1 1])
if noiseon
    plot(signal.t, signal.snom, '-','LineWidth',1,'Color',[0 0 0 0.6]);
    l=legend('$x(t)$','$x(kT)$','$\hat{x}(t)$','$\tilde{x}(t)$','$\bar{x}(t)$');
else
    l=legend('$x(t)$','$x(kT)$','$\hat{x}(t)$','$\tilde{x}(t)$');
end
set(l,'interpreter','latex');
xlabel('Time (s)');
ylabel('Signal Value')

figure(11111);hold on;
plot(signal.t, signal.s, '-','LineWidth',4,'Color',[0.8 0.8 0.8 0.4]);
stem(sampled.t3, sampled.s3, '-o','LineWidth',1,'Color',[0 1 1 0.5]);
plot(sampled.t3, sampled.s3,'-','LineWidth',3,'Color',[0 1 1 1])
plot(signal.t, sampled.p3,'--','LineWidth',1,'Color',[0 1 1 1])
if noiseon
    plot(signal.t, signal.snom, '-','LineWidth',1,'Color',[0 0 0 0.6]);
    l=legend('$x(t)$','$x(kT)$','$\hat{x}(t)$','$\tilde{x}(t)$','$\bar{x}(t)$');
else
    l=legend('$x(t)$','$x(kT)$','$\hat{x}(t)$','$\tilde{x}(t)$');
end
set(l,'interpreter','latex');
xlabel('Time (s)');
ylabel('Signal Value')


%% ----------------- Simple Sampling Example ------------------------------
clear all
close all
clc

signalmaxtime = 0.1;
signalfreq    = 25;
signaltime    = 0:1/1000:signalmaxtime;
signalsimp    = @(f,t) cos(2*pi*f*t);
samplefreq(1) = (4/5)*signalfreq;   % < Nyq Sampling Frequencies
samplefreq(2) = 2*signalfreq;       % @ Nyq Sampling Frequencies
samplefreq(3) = 2*(2*signalfreq);   % > Nyq Sampling Frequencies

% --------------------- Low Sampling Frequencies --------------------------
for k=1:1:length(samplefreq)
    fs = samplefreq(k);
    fa = (abs(signalfreq-fs));
    fb = (abs(signalfreq+fs));

    % Sample Signal
    sampletime = 0:1/fs:signalmaxtime;          % Sample Times (vector)          

    a1 = signalsimp(signalfreq,sampletime);     % Sample Values (Freq A)  
    a2 = signalsimp(signalfreq,sampletime);     % Sample Values (vector)

    s1 = signalsimp(fa,signaltime);             % Fitted Sample Values (Freq A)
    s2 = signalsimp(fb,signaltime);             % Fitted Sample Values (Freq B) 

    % Show Aliasing and Sampling
    figure(2+k);hold on;
    plot(signaltime, signalsimp(signalfreq,signaltime), ...
        '-','LineWidth',4,'Color',[0.8 0.8 0.8]);
    plot(signaltime, s1, 'k-','LineWidth',1);
    plot(signaltime, s2, 'k--','LineWidth',1);  
    plot(sampletime, a1, 'ko','MarkerSize',10);
    plot(sampletime, a2, 'k.','MarkerSize',10);
    l=legend('$x(t)$','$\hat{x}_a(t)$','$\hat{x}_b(t)$','$x_a(kT)$','$x_b(kT)$');
    set(l,'interpreter','latex');
    xlabel('Time (s)');
    ylabel('Signal Value')
end


