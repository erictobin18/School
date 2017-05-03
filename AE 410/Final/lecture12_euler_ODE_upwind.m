function dqdt = euler_ODE_upwind(t,q,Nx,gamma,D_plus,D_minus);

rho = q(1     :  Nx);
u   = q(Nx+1  :2*Nx); 
p   = q(2*Nx+1:3*Nx);

d_rho_plus  = D_plus*rho;
d_rho_minus = D_minus*rho;

d_u_plus    = D_plus*u;
d_u_minus   = D_minus*u;

d_p_plus    = D_plus*p;
d_p_minus   = D_minus*p;

for m = 1:Nx;
  c(m) = sqrt(gamma*p(m)/rho(m));

  X = [ 1    1/c(m)^2         1/c(m)^2          ;
        0   -1/(rho(m)*c(m))  1/(rho(m)*c(m))   ;
        0    1                1                 ];  

  Lambda = diag([u(m) u(m)-c(m) u(m)+c(m)]);    

  Lambda_plus  =  subplus( Lambda);
  Lambda_minus = -subplus(-Lambda);

  d_w_plus  = inv(X)*[d_rho_plus(m);  d_u_plus(m);  d_p_plus(m)];
  d_w_minus = inv(X)*[d_rho_minus(m); d_u_minus(m); d_p_minus(m)];
      
  d_w  = Lambda_plus*d_w_plus + Lambda_minus*d_w_minus;
  d_q  = X*d_w;

  dqdt(m     ,1)  = - d_q(1);
  dqdt(m+Nx  ,1)  = - d_q(2);
  dqdt(m+2*Nx,1)  = - d_q(3);
end
