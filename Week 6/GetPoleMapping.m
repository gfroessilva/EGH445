%% EGH445 - Modern Control
%--------------------------------------------------------------------------
% Modern Control
% **** Support Functions
%     **** Pole Mapping
%
% Author: Aaron McFadyen
%--------------------------------------------------------------------------

function [mapping] = GetPoleMapping(A,B,C,D,T)


% Find Open Loop Poles (s-z mapping)
css = ss(A,B,C,D);                       % Continuous State Space Models
[num,den] = ss2tf(A,B,C,D);
Gs.tf = tf(num,den);                    % Continuous Transfer Function
Gs.poles = roots(den);                  % Continuous System Poles

for k=1:1:length(T)
    % Recall z = e^(Ts)
    poles(k,:) = exp(T(k).*ones(size(Gs.poles)).*Gs.poles)'; % s -> z mapping
    disp(horzcat('Discrete Sampling Time/Period: ',num2str(T(k)),' seconds'));
    for i=1:1:length(poles(k,:))
        disp(horzcat('Continuous (s) Pole ',num2str(i),': ',num2str(Gs.poles(i,1))));
        disp(horzcat('Discrete   (z) Pole ',num2str(i),': ',num2str(poles(k,i))));
    end
    disp('')
end
mapping.contpoles = Gs.poles;
mapping.discpoles = poles;
mapping.period = T;

end

