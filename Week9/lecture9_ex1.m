G = [0.5 0.1; 0 0.5]; H = [1 0; 0 1];
syms k1 k2 k3 k4 real 
K = [k1 k2; k3 k4];

syms z

Gcl = G - H*K;
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