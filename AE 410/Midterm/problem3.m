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
  N = 10;%500;
  L  = 1;
  
  % time interval 
  t_max = 0.05;
  N_t   = 100;
  dt    = t_max/N_t;
  t_plot = [0:dt:t_max]; 

% --------------------------------------------

% construct mesh
dx = L/N;
x  = dx*linspace(0,N-1,N).';

% construct FD matrix
PeriodicFlag = 1;
Dp = sparse(upwind_3rd_order(x,PeriodicFlag));
Dm = sparse(downwind_3rd_order(x,PeriodicFlag));
% initial condition
rho(1:N,1) = 1.0;
  u(1:N,1) = 0.0;
  p(1:N,1) = 0.1 + exp(-((x-0.5)/0.05).^2);

E = p./((gamma-1).*rho) + 1/2*(u.^2);

q_0(1:N,1)       = rho;
q_0(N+1:2*N,1)   = u;
q_0(2*N+1:3*N,1) = p;

[t,q] = ode45(@(t,q) problem3_system(t,q,N,gamma,Dp,Dm),t_plot,q_0);

figure;
for n = 1:length(t)
  rho   = q(n,1      :N);
  u = q(n,N+1  :2*N); 
  p = q(n,2*N+1:3*N);

  E = p./((gamma-1).*rho) + 1/2*(u.^2);
  e   = E - 1/2*(u.^2);        

  subplot(2,2,1); plot(x,rho); ylabel('\rho');  ylim([0 5]);
  subplot(2,2,2); plot(x,u);   ylabel('u');     ylim([-1 1]); 
  subplot(2,2,3); plot(x,p);   ylabel('p');     ylim([0 1]);
  subplot(2,2,4); plot(x,e);   ylabel('e');     ylim([0 3]);
  drawnow;
end
