function [ grid ] = Setup_Quadrants( x, y, p, rho, u, v, gamma )
nx = length(x);
ny = length(y);
n2xc = ceil(nx/2);
n2xf = floor(nx/2);
n2yc = ceil(ny/2);
n2yf = floor(ny/2);

s1 = ones(n2xc,n2yc);
s2 = ones(n2xc,n2yf);
s3 = ones(n2xf,n2yc);
s4 = ones(n2xf,n2yf);

rho_E = ideal_rho_E(p(:),rho(:),rho(:).*u(:),rho(:).*v(:),gamma);
rho_E = reshape(rho_E,2,2);

grid = zeros(length(x),length(y),4);
grid(:,:,1) = [
    rho(1,1)*s1 rho(1,2)*s3;
    rho(2,1)*s2 rho(2,2)*s4;
    ];
grid(:,:,2) = [
    rho(1,1)*u(1,1)*s1 rho(1,2)*u(1,2)*s3;
    rho(2,1)*u(2,1)*s2 rho(2,2)*u(2,2)*s4;
    ];
grid(:,:,3) = [
    rho(1,1)*v(1,1)*s1 rho(1,2)*v(1,2)*s3;
    rho(2,1)*v(2,1)*s2 rho(2,2)*v(2,2)*s4;
    ];
grid(:,:,4) = [
    rho_E(1,1)*s1 rho_E(1,2)*s3;
    rho_E(2,1)*s2 rho_E(2,2)*s4;
    ];

grid = permute(grid,[3,1,2]);
end

function out = ideal_rho_E(p,rho,rho_u,rho_v,gamma) % solve equation of state for rho*E
rho_norm_v_squared = rho_u.^2 + rho_v.^2;
out = (p/(gamma - 1)) + (1/2)*(rho_norm_v_squared.^2)./rho;
end