function [dx] = XYmass(x,T,u)
A = [1 T 0 0  
     0 1 0 0  
     0 0 1 T
     0 0 0 1];
 B = [0 0     
      1 0
      0 0
      0 1];
 dx = zeros(4,1);
 dx = A*x+T*B*u;
end
