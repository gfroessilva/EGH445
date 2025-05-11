function [x,P] = Fkalman(x,P,A_d,B_d,C_d,T,R,Q,z,w)
%prediçao
xant= A_d*x+B_d*w;
Pk= A_d*P*A_d' + R;
%correçao
K=Pk*C_d' * inv(C_d*Pk*C_d'+Q);
x=xant+K*(z-C_d*xant);
P=(eye(4)-K*C_d)*Pk;
     
end

