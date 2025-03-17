%% EGH445 - Modern Control
%--------------------------------------------------------------------------
% Modern Control
% **** Support Functions
%     **** System Models
%
% Author: Aaron McFadyen
%--------------------------------------------------------------------------


%% Plant Selection Function 

function [A,B,C,D] = GetPlantModel(type)

switch type
    
    case 'mass-spring'    % Example Lecture Slides 
        % Continuous State-Space    
        A = [0 1;-1 -2];                    % Plant Matrix  (C/s)  
        B = [0;1];                          % Input Matrix  (C/s)
        C = [0 1];                          % Output Matrix (C/s)
        %C = [1 0];
        D = 0;                              % Output Control Matrix (C/s) [use 0 for all mtx sizes]
    
    case 'nice'    % Example 3-3 Nice (p.138)
        % Continuous Time Transfer Function
        num         = [0 2 1];              % Transfer Function Numerator           (C/s)
        den         = [1 7 9];              % Transfer Function Denominator         (C/s)
        [A,B,C,D]  	= tf2ss(num,den);     	% Continuous State-Space (Controller Can. Form)
    

end



