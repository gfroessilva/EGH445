function jacob = calc_jacob(f, x0, dx)
%% Help for calc_jacob
%    Input syntax:
%    -- f:  N-dimentional functions with respect to x
%    -- x0: N-dimentional vector for x for Jacobian calc
%    -- dx: vector length

%% Coding Area
% we reshape x0 to vertical vector
x0 = x0(:);

% Calc the Jacobian
dim = numel(x0);
jacob = [];

for cnt = 1:dim
    % we create a vector x0+dx
   partial_x = zeros(size(x0));
   partial_x(cnt) = dx;
   
   % calc partial f[1...n] to [0...dx_i...0]
   partial_f =  f(x0+partial_x) - f(x0);
   
   % reshape partial f to a vector
   partial_f = partial_f(:);
   
   % assign jacobian column
   jacob(:,cnt) = partial_f/dx;
end

end
