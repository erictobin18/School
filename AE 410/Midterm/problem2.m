clear all; clc;

% -------------------------------------------
  % include path for matrix operator m files
  addpath ./matrix_operators

  % set default plot attributes
  default_plot_attributes  

  % spatial discretization
  N = 201;
  L  = 2;
  dx = L/(N-1);
  x  = dx*linspace(-(N-1)/2,(N-1)/2,N).';

  % time interval 
  t_max = 2;
  N_t   = 200;
  dt    = t_max/N_t;
  t_plot = [0:dt:t_max];

  % initial condition
  u = [ones((N-1)/2,1)*1;ones((N+1)/2,1)*.5];
  v = [ones((N-1)/2,1)*.5;ones((N+1)/2,1)*.7];	
  w = [ones((N-1)/2,1)*1;ones((N+1)/2,1)*.2];
  
  % Viscosity
  viscosity = dx/10;
% --------------------------------------------

% construct FD matrix
PeriodicFlag = 1;
D = central_2nd_order(x,PeriodicFlag);
D2 = viscosity*central_2nd_order_second_derivative(x,PeriodicFlag);
Dp = upwind_3rd_order(x,PeriodicFlag);
Dm = downwind_3rd_order(x,PeriodicFlag);

q_0 = [u; v; w];
   
% Artificial viscosity solution using ode45
[t1,q1] = ode45(@(t1,q1) problem2_system(t1,q1,D,D2),t_plot,q_0);

Lp = diag([(1+sqrt(5))/2;1;0]);
Lm = diag([0;0;(1-sqrt(5))/2]);
X = [1/3,1,1/3;-2/(3*(1-sqrt(5))),1,-2/(3*(1+sqrt(5)));1,1,1];
Ap = inv(X)*(Lp*X);
Am = inv(X)*(Lm*X);

% Upwind solution using ode45
[t2,q2] = ode45(@(t2,q2) problem2_system_upwind(t2,q2,Ap,Am,Dp,Dm),t_plot,q_0);
clf;

figure(1);
contour(q1(:,1:N),'k');
xlabel('x');
ylabel('t');
title('Viscosity Rho');

figure(2);
contour(q1(:,N+1:2*N),'k');
xlabel('x');
ylabel('t');
title('Viscosity Velocity');

figure(3);
contour(q1(:,2*N+1:end),'k');
xlabel('x');
ylabel('t');
title('Viscosity Pressure');

figure(4);
contour(q2(:,1:N),'k');
xlabel('x');
ylabel('t');
title('Upwind Rho');

figure(5);
contour(q2(:,N+1:2*N),'k');
xlabel('x');
ylabel('t');
title('Upwind Velocity');

figure(6);
contour(q2(:,2*N+1:end),'k');
xlabel('x');
ylabel('t');
title('Upwind Pressure');

figure
for n = 1:N_t	
	t = n*dt;

  u = q2(n,1:N);
  v = q2(n,N+1:2*N);
  w = q2(n,2*N+1:3*N);

  subplot(1,3,1); 
    plot(x,u); ylim([-1 1]); xlabel('x'); ylabel('u_1');
  subplot(1,3,2); 
    plot(x,v); ylim([-1 1]); xlabel('x'); ylabel('u_2');
  subplot(1,3,3); 
    plot(x,w); ylim([-1 1]); xlabel('x'); ylabel('u_3');
  
  pause(0.001);

	drawnow;
end 
