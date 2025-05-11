
set(0, "DefaultAxesFontSize", 18)
set(0, "DefaultAxesLineWidth", 3)
%% System
A = [0 1; -1 -2]; B = [0; 1]; C = [1 0]; D = 0;

%% Continuous Controller
K_c = 2.8;
Gc_s = tf(K_c*[1 2], [1 4]);
G.Numerator = Gc_s.Numerator{1};
G.Denominator = Gc_s.Denominator{1};

sysC = ss(A,B,C,D);
t = 0:0.1:10;
[y,t] = lsim(feedback(series(sysC, Gc_s),1), ones(size(t)), t);

% calculate input u(t)
e = 1 - y; % Calculate the error signal
u = lsim(Gc_s, e, t);

f = figure(1); clf; f.Theme = 'light';
subplot(2,1,1); plot(t, y, linewidth=4); grid on; title('Output y(t)', fontsize=20);
hold on
axis([0 10 0 1])
xlabel('Time (s)'); ylabel('y(t)');
legend('Continous Output')
subplot(2,1,2); plot(t, u, linewidth=4); grid on; title('Input u(t)', fontsize=20);
hold on
xlabel('Time (s)'); ylabel('u(t)');
axis([0 10 0 3])
legend('Continous Input')

T = 1;
%% Tustin
Gd_tustin = c2d(Gc_s, T, 'tustin');
disp('Gd_tustin (T=1.0s):');


Cd = Gd_tustin;

out = sim('EmulationDesign');

% plot results

subplot(2,1,1); 
plot(out.y, linewidth=4, DisplayName="Tustin Output"); 
grid on; title("Output y(t), T="+T, fontsize=20);
axis([0 10 0 1])
xlabel('Time (s)'); ylabel('y(t)');
subplot(2,1,2); 
plot(out.u, linewidth=4, DisplayName="Tustin Input"); 
grid on; title('Input u(t)', fontsize=20);
xlabel('Time (s)'); ylabel('u(t)');
axis([0 10 0 3])

%% ZOH

K_c = 2.8;
Gc_s = tf(K_c*[1 2], [1 4]);
Gd_zoh = c2d(Gc_s, T, 'zoh');

Cd = Gd_zoh;

out = sim('EmulationDesign');

% plot results

subplot(2,1,1); 
plot(out.y, linewidth=4, DisplayName="ZOH Output"); 
grid on; title("Output y(t), T="+T, fontsize=20);
axis([0 10 0 1])
xlabel('Time (s)'); ylabel('y(t)');
subplot(2,1,2); 
plot(out.u, linewidth=4, DisplayName="ZOH Input"); 
grid on; title('Input u(t)', fontsize=20);
xlabel('Time (s)'); ylabel('u(t)');
axis([0 10 0 3])

%% MPZ

K_c = 2.8;
Gc_s = tf(K_c*[1 2], [1 4]);
Gd_zoh = c2d(Gc_s, T, 'mpz');

Cd = Gd_zoh;

out = sim('EmulationDesign');

% plot results

subplot(2,1,1); 
plot(out.y, linewidth=4, DisplayName="ZOH Output"); 
grid on; title("Output y(t), T="+T, fontsize=20);
axis([0 10 0 1])
xlabel('Time (s)'); ylabel('y(t)');
subplot(2,1,2); 
plot(out.u, linewidth=4, DisplayName="ZOH Input"); 
grid on; title('Input u(t)', fontsize=20);
xlabel('Time (s)'); ylabel('u(t)');
axis([0 10 0 3])