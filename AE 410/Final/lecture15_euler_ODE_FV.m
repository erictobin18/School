function dqdt = euler_ODE_FV(t,q,N,dx,gamma);
	
	rho   = q(1    :  N);
	rho_u = q(N+1  :2*N);	
	rho_E = q(2*N+1:3*N);
	
	u = rho_u./rho;
	p = (gamma-1)*(rho_E - 1/2*rho.*u.^2);
	
	% store temporary variables 
	f1 =  rho_u;
	f2 =  rho_u.*u + p;
	f3 = (rho_E    + p).*u;


  % Compute F_{i+1/2} = 0.5*(F_{i+1/2}^L + F_{i+1/2}^R) -----------------------
	for i = 1:N-1
		f1_plus_half(i) = 0.5*(f1(i) + f1(i+1));
		f2_plus_half(i) = 0.5*(f2(i) + f2(i+1));
		f3_plus_half(i) = 0.5*(f3(i) + f3(i+1));
	end
	
	% assuming periodic BCs, i.e. F_{N+1/2}^R = F_{1-1/2}^R
	i = N;	
	f1_plus_half(i) = 0.5*(f1(i) + f1(1));
	f2_plus_half(i) = 0.5*(f2(i) + f2(1));
	f3_plus_half(i) = 0.5*(f3(i) + f3(1));
	% ---------------------------------------------------------------------------


	% Compute F_{i-1/2} = 0.5*(F_{i-1/2}^L + F_{i-1/2}^R) -----------------------
	for i = 2:N
		f1_minus_half(i) = 0.5*(f1(i-1) + f1(i));
		f2_minus_half(i) = 0.5*(f2(i-1) + f2(i));
		f3_minus_half(i) = 0.5*(f3(i-1) + f3(i));
	end
	
	% assuming periodic BCs
	i = 1;	
	f1_minus_half(i) = 0.5*(f1(N) + f1(i));
	f2_minus_half(i) = 0.5*(f2(N) + f2(i));
	f3_minus_half(i) = 0.5*(f3(N) + f3(i));
	% ---------------------------------------------------------------------------

  dqdt = zeros(3*N,1);
	for i = 1:N;	
		dqdt(i    ,1) = - (f1_plus_half(i) - f1_minus_half(i))/dx;
		dqdt(N+i  ,1) = - (f2_plus_half(i) - f2_minus_half(i))/dx;
		dqdt(2*N+i,1) = - (f3_plus_half(i) - f3_minus_half(i))/dx;
	end
