function out = ideal_rho_E(p,rho,rho_u,gamma) % solve equation of state for rho*E
out = (p/(gamma - 1)) + (1/2)*(rho_u.^2)./rho;
end