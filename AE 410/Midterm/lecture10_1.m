clear all; clc;

% -------------------------------------------
  % include path for matrix operator m files
  addpath ./matrix_operators

  % set default plot attributes
  default_plot_attributes  

  % spatial discretization
  N = 101;
  L  = 1;
  dx = L/N;
  x  = dx*linspace(0,N-1,N).';

  % time interval 
  t_max = 2;
  N_t   = 200;
  dt    = t_max/N_t;
  t_plot = [0:dt:t_max];

  % initial condition
  u_1 = zeros(N,1);
  u_2 = exp(-((x-0.5)/0.1).^2);	
  u_3 = zeros(N,1);
% --------------------------------------------

% construct FD matrix
PeriodicFlag = 1;
D = central_2nd_order(x,PeriodicFlag);

q_0 = [u_1; u_2; u_3];
   
% Solution using ode45
[t,q] = ode45(@(t,q) lecture10_1_hyperbolic_system_of_pdes(t,q,D),t_plot,q_0);

figure
for n = 1:N_t	
	t = n*dt;

  u_1 = q(n,1:N);
  u_2 = q(n,N+1:2*N);
  u_3 = q(n,2*N+1:3*N);

  subplot(1,3,1); 
    plot(x,u_1); ylim([-1 1]); xlabel('x'); ylabel('u_1');
  subplot(1,3,2); 
    plot(x,u_2); ylim([-1 1]); xlabel('x'); ylabel('u_2');
  subplot(1,3,3); 
    plot(x,u_3); ylim([-1 1]); xlabel('x'); ylabel('u_3');
  
  pause(0.001);

	drawnow;
end 
