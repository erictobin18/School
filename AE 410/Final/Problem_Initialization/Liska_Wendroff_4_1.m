function [ out ] = Liska_Wendroff_4_1( x, y, gamma )
% Case 3

p = [.3,1.5;.029,.3];
rho = [.5323,1.5;.138,.5323];
u = [1.206,0;1.206 0];
v = [0,0;1.206,1.206];

out = Setup_Quadrants(x, y, p, rho, u, v, gamma);
end

