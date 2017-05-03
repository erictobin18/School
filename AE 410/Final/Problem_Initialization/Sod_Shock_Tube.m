function [ out ] = Sod_Shock_Tube( x, gamma )
% the three columns are columns for rho, rho*u, rho*E
n = length(x);
n2 = floor(n/2);
rho1 = 1;
p1 = 1;
u1 = 0;
rho2 = .125;
p2 = .1;
u2 = 0;
out = [ 
    ones(n2,1)*rho1   ones(n2,1)*rho1*u1   ones(n2,1)*ideal_rho_E(p1,rho1,rho1*u1,gamma);
    ones(n-n2,1)*rho2 ones(n-n2,1)*rho2*u2 ones(n-n2,1)*ideal_rho_E(p2,rho2,rho2*u2,gamma);
    ]';
end
