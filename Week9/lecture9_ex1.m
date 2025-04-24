clear
G = [1 0.1; 0 0.8]; H = [1 0; 0 1];

z1 = 0.4;
z2 = 0.6;

syms z
syms k_1 k_2 k_3 k_4 real
K = [k_1 k_2; k_3 k_4];

desiredPoly = collect((z-z1)*(z-z2));
%%


Gcl = G - H*K;

z*eye(2) - Gcl;

actualPoly = vpa(collect(det(z*eye(2) - Gcl)),3);

eq1 = (k_1 + k_4 - 1.8) == -1;
eq2 = -0.8*k_1 + 0.1*k_3 - 1.0*k_4 + k_1*k_4 - 1.0*k_2*k_3 + 0.8 == 0.24;
%%
k2 = 0; k3 = 0;
eqq2 = subs(eq2, [k_3, k_2], [k3, k2]);
sol = solve([eq1, eqq2], [k_1 k_4]);
k1 = double(sol.k_1);
k4 = double(sol.k_4);

K1 = [k1(1) k2; k3 k4(1)]
eigGcl = eig(G-H*K1)
K2 = [k1(2) k2; k3 k4(2)]
eigGcl = eig(G-H*K2)


%%
charEq = det(z*eye(2) - Gcl);

eq1 = k1 + k4 - 1.0 == -1.1;
eq2 = k1*k4 - 0.5*(k1 + k4) + 0.1*k3 - k2*k3 +0.25 == 0.3;

%%
home
disp('Control Gain K1:')
k2_ = 0; k3_ = 0;
% k1 = -0.29          1.2
% k4 = -0.19 1.3

eq11 = subs(eq1, [k2 k3], [k2_ k3_]);
eq21 = subs(eq2, [k2 k3], [k2_ k3_]);
K = solve([eq11; eq21], [k1 k4]);
k1_ = vpa(K.k1,2)
k2_
k3_
k4_ = vpa(K.k4,2)

K1 = [k1_(1) k2_; k3_ k4_(1)];
K2 = [k1_(2) k2_; k3_ k4_(2)];
%%
disp('Control Gain K2:')
k2_ = 1; k3_ = 1.25;
eq11 = subs(eq1, [k2 k3], [k2_ k3_]);
eq21 = subs(eq2, [k2 k3], [k2_ k3_]);
K = solve([eq11 eq21], [k1 k4]);
k1_ = vpa(K.k1,2)
k2_
k3_ 
k4_ = vpa(K.k4,2)

K2 = [k1_ k2_; k3_ k4_];