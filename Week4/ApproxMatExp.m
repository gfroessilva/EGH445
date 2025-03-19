%% EGH445 - Modern Control
%--------------------------------------------------------------------------
% Modern Control
% **** Support Functions
%     **** Matrix Exp Approximation
%
% Author: Aaron McFadyen
%--------------------------------------------------------------------------


%% Function Definition

function [G, H] = ApproxMatExp(A, B, T, n, showsolution)
% A = Continuous time plant matrix (A)
% A = Continuous time input matrix (B)
% T = Sampling Time
% n = Order of approximation for matrix exponential

G = eye(size(A));
H = zeros(size(B));

for i=1:1:n
    G = G + A^i.*T^i.*(1/factorial(i));
    %------ STUDENTS: Insert code to approximate H ------
    % H = ?
    %----------------------------------------------------
    if showsolution
        disp(horzcat('Matrix Exponetial Approximation - Order ',num2str(i),':'))
        disp('G: '); G
        disp('H: '); H
    end
end
