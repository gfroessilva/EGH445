function [u] = controle(x,r)
    kp = .05;
    kd =-.5;
    u = [(r(1)-x(1))*kp + x(2)*kd; (r(2)-x(3))*kp + x(4)*kd];
end
