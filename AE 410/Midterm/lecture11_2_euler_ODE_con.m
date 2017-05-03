function dqdt = lecture11_2_euler_ODE_con(t,q,N,gamma,D);

rho   = q(1    :  N);
rho_u = q(N+1  :2*N); 
rho_E = q(2*N+1:3*N);

u = rho_u./rho;
p = (gamma-1)*(rho_E - 1/2*rho.*u.^2);

dqdt(1    :  N,1) = - D*(rho_u          );
dqdt(N+1  :2*N,1) = - D*(rho_u.*u + p   );
dqdt(2*N+1:3*N,1) = - D*(rho_E.*u + p.*u);

