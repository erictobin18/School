clear all; clc;

% -------------------------------------------
  % include path for matrix operator m files
  addpath ./matrix_operators

  % set default plot attributes
  default_plot_attributes  

  % wave speed
  a = 1;

  % spatial discretization
  N = 201;
  L  = 1;
  dx = L/N;
  x  = dx*linspace(0,N-1,N).';

  % time interval 
  t_max = 2;
  N_t   = 200;
  dt    = t_max/N_t;
  t_plot = [0:dt:t_max];

  % initial condition
  u_1_0 = zeros(N,1);
  u_2_0 = exp(-((x-0.5)/0.1).^2);	
  u_3_0 = zeros(N,1);
% --------------------------------------------

% construct FD matrix
PeriodicFlag = 1;
D = central_2nd_order(x,PeriodicFlag);

A = [2  4  0; 
     0 -2  1; 
     0  1  2];

% compute eigenvalues and eigenvectors of A
[X,Lambda] = eig(A);
lambda = diag(Lambda);  

% compute inverse of X
X_inv = inv(X);

% initial condition in characteristic variables,  w = X^{-1}*u
w_1_0 = X_inv(1,1)*u_1_0 + X_inv(1,2)*u_2_0 + X_inv(1,3)*u_3_0;
w_2_0 = X_inv(2,1)*u_1_0 + X_inv(2,2)*u_2_0 + X_inv(2,3)*u_3_0;
w_3_0 = X_inv(3,1)*u_1_0 + X_inv(3,2)*u_2_0 + X_inv(3,3)*u_3_0;

% Solution using ODE45
[t,w_1] = ode45(@(t,w_1) lecture10_2_diagonalized(t,w_1,-lambda(1)*D),t_plot,w_1_0);
[t,w_2] = ode45(@(t,w_2) lecture10_2_diagonalized(t,w_2,-lambda(2)*D),t_plot,w_2_0);
[t,w_3] = ode45(@(t,w_3) lecture10_2_diagonalized(t,w_3,-lambda(3)*D),t_plot,w_3_0);

% solution in physical coordinates, u = X*w
u_1 = X(1,1)*w_1 + X(1,2)*w_2 + X(1,3)*w_3;
u_2 = X(2,1)*w_1 + X(2,2)*w_2 + X(2,3)*w_3;
u_3 = X(3,1)*w_1 + X(3,2)*w_2 + X(3,3)*w_3;

figure
for n = 1:N_t	
	t = n*dt;

  % characteristic variables
  subplot(2,3,1); 
    plot(x,w_1(n,:)); ylim([-1.5 1.5]); xlabel('x'); ylabel('w_1');
  subplot(2,3,2); 
    plot(x,w_2(n,:)); ylim([-1.5 1.5]); xlabel('x'); ylabel('w_2');
  subplot(2,3,3); 
    plot(x,w_3(n,:)); ylim([-1.5 1.5]); xlabel('x'); ylabel('w_3');  

  % physical variables
  subplot(2,3,4); 
    plot(x,u_1(n,:)); ylim([-1.5 1.5]); xlabel('x'); ylabel('u_1');
  subplot(2,3,5); 
    plot(x,u_2(n,:)); ylim([-1.5 1.5]); xlabel('x'); ylabel('u_2');
  subplot(2,3,6); 
    plot(x,u_3(n,:)); ylim([-1.5 1.5]); xlabel('x'); ylabel('u_3');
  
  pause(0.1);

	drawnow;
end 
