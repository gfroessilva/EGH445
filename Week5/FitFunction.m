%% EGH445 - Modern Control
%--------------------------------------------------------------------------
% Modern Control
% **** Support Functions
%     **** Signal Reconstruction from Samples
%
% Author: Aaron McFadyen
%--------------------------------------------------------------------------


%% Function Definition
function [tp, yp] = FitFunction(t, y)

yu = max(y);
yl = min(y);
yr = (yu-yl);                               % Range of ‘y’
yz = y-yu+(yr/2);
zx = t(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
per = 2*mean(diff(zx));                     % Estimate period
ym = mean(y);                               % Estimate offset

fit = @(p,t)  p(1).*cos(2*pi*t.*p(2)) +  p(3).*sin(2*pi*t.*p(4));  	% Function to fit
fcn = @(p) sum((fit(p,t) - y).^2);                                  % Least-Squares cost function
s = fminsearch(fcn, [1 1/per 1 1/per]) ;                          	% Minimise Least-Squares
%s = fminsearch(fcn, [yr;  per;  -1;  ym]) ;

tp = linspace(min(t),max(t),1000);
yp = fit(s,tp);

% figure(222); hold on
% plot(t,y,'b',  tp, yp, 'r');
% grid;



end

