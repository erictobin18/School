clear all;
clc;

%=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1D Euler equations
%
% dq/dt + dF(q)/dx = 0
%
% where
%
%   q = [q_1] = [rho  ],  F(Q) = [Q_2       ]
%       [q_2]   [rho*u]          [Q_2*u + p ]
%       [q_3]   [rho*E]          [(Q_3*u+p)u]
%
% With the state equation:
%
%   p = (\gamma -1)(\rho*E - 1/2*\rho*u^2)  
%
%=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------
  % include path for matrix operator m files
  addpath ./matrix_operators

  % set default plot attributes
  default_plot_attributes  

  % ratio of specific heats 
  gamma = 1.4;

  % discretization parameters
  N = 500;
  L  = 1;
  
  % time interval 
  t_max = 0.3;
  N_t   = 100;
  dt    = t_max/N_t;
  t_plot = [0:dt:t_max]; 

% --------------------------------------------

% construct mesh
dx = L/N;
x  = dx*linspace(0,N-1,N).';

% initial condition
rho(1:N,1) = 1.0;
  u(1:N,1) = 0.0;
  p(1:N,1) = 0.1 + exp(-((x-0.5)/0.05).^2);

E = p./((gamma-1).*rho) + 1/2*(u.^2);

q_0(1:N,1)       = rho;
q_0(N+1:2*N,1)   = rho.*u;
q_0(2*N+1:3*N,1) = rho.*E;

[t,q] = ode45(@(t,q) lecture15_euler_ODE_FV(t,q,N,dx,gamma),t_plot,q_0);

figure;
for n = 1:length(t)
  rho   = q(n,1      :N);
  rho_u = q(n,N+1  :2*N); 
  rho_E = q(n,2*N+1:3*N);

  u   = rho_u./rho;
  E   = rho_E./rho;
  e   = E - 1/2*(u.^2); 
  p   = e.*((gamma-1)*rho);       

  subplot(2,2,1); plot(x,rho); ylabel('\rho');  ylim([0 5]);
  subplot(2,2,2); plot(x,u);   ylabel('u');     ylim([-1 1]); 
  subplot(2,2,3); plot(x,p);   ylabel('p');     ylim([0 1]);
  subplot(2,2,4); plot(x,e);   ylabel('e');     ylim([0 3]);
  drawnow;
end
