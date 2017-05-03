function [ out ] = Liska_Wendroff_4_2( x, y, gamma )
% Case 12

p = [1.0,.4;1.0,1.0];
rho = [1.0,.5313;.8,1.0];
u = [.7276,0;0,0];
v = [0,0;0,.7276];

out = Setup_Quadrants(x, y, p, rho, u, v, gamma);
end

