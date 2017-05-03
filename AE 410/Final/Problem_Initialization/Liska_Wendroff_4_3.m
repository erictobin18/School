function [ out ] = Liska_Wendroff_4_3( x, y, gamma )
% Case 15

p = [.4, 1.0;.4, .4];
rho = [.5197,1.0;.8,.5313];
u = [-.6259,.1;.1,.1];
v = [-.3,-.3;-.3,.4276];

out = Setup_Quadrants(x, y, p, rho, u, v, gamma);
end

