A = [0 1; 0 0]; B = [0;1]; C = eye(2); D = 0;

sysc = ss(A,B,C,D);

Ts = 0.1;
sysd = c2d(sysc, Ts, 'zoh');

[G,H,C,D] = ssdata(sysd);

Q = [30 0; 0 10]; R = 5;

[K, P, E] = dlqr(G,H,Q,R);

t = 0:0.1:10;

sysd_cl = ss(G-H*K, H, C, D, Ts);

y = initial(sysd_cl, [1; 0], t);

x1 = y(:,1);
x2 = y(:,2);

stairs(t, x1, 'b', 'linewidth', 3); hold on
stairs(t, x2, 'g', 'linewidth', 3);


u = -K*y';
stairs(t, u, 'r', 'LineWidth', 3);
legend('x1', 'x2', 'u', fontsize=16)
a = gca; a.FontSize=16;
grid on
hold off