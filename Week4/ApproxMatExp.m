%% EGH445 - Modern Control
%--------------------------------------------------------------------------
% Digital Control
% ****Support Functions
%     ****Matrix Exponential Approximation (nth order)
%
% Author:   Dr Aaron McFadyen
% Created:  August 2018
% Modified: August 2018
%--------------------------------------------------------------------------


%% Function Definition

function [G] = ApproxMatExp(A, T, n, showsolution)
% A = Continuous time plant matrix (A)
% T = Sampling Time
% n = Order of approximation for matrix exponential

G = eye(size(A));

for i=1:1:n
    G = G + A^i.*T^i.*(1/factorial(i)); 
    if showsolution
        disp(horzcat('Matrix Exponetial Approximation - Order ',num2str(i),':'))
        G
    end
end


% Resources
%G3 = eye(size(S.A)) + S.A.*T + S.A^2.*T^2.*(1/factorial(2)) + ...
%    S.A^3.*T^3.*(1/factorial(3))