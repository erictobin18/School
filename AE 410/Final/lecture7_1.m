clear all; clc;

%=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Linear convection equation with with periodic boundary conditions
%
%	  du/dt + du/dx = 0,  x \in [0,1),  u(x,0) = u_0(x)
%
% 3rd order upwind and downwind FD; see Eq 7.1a and Eq 7.1b in lecture notes.
%
%=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------
	N   = 51;
	L   = 1;
	dx  = L/N;
	x   = dx*linspace(0,N-1,N).';
	
	% default plot attributes
	set(0,'defaultfigurecolor',		  [1 1 1]);
	%set(0,'defaultfigureposition',	[10 700 800 450])
% --------------------------------------------

% central, 2nd-order accurate
stencil = [-1 0 1];
for i = 2:N-1
  D(i,i-1:i+1) = stencil;
end
D(1,  1:2) = stencil(2:3); D(1,N) = stencil(1);
D(N,N-1:N) = stencil(1:2); D(N,1) = stencil(3);

D_central_2nd = D./(2*dx); clear D;

% upwind, 3rd-order accurate
stencil = [1 -6 3 2];
for i = 3:N-1
  D(i,i-2:i+1) = stencil;
end
D(1,  1:2) = stencil(3:4);   D(1,N-1:N) = stencil(1:2);
D(2,  1:3) = stencil(2:4);   D(2,    N) = stencil(1);
D(N,N-2:N) = stencil(1:3);   D(N,    1) = stencil(4);

D_upwind_3rd = D./(6*dx); clear D;

% downwind, 3rd-order accurate
stencil = [-2 -3 6 -1];
for i = 2:N-2
  D(i,i-1:i+2) = stencil;
end
D(1,    1:3) = stencil(2:4);   D(1,    N) = stencil(1);
D(N-1,N-2:N) = stencil(1:3);   D(N-1,  1) = stencil(4);
D(N,  N-1:N) = stencil(1:2);   D(N  ,1:2) = stencil(3:4);

D_downwind_3rd = D./(6*dx); clear D;

% choose D
D = D_central_2nd;

% analytical eigenvalues/vectors (Eq. 6.10 and 6.12 in lecture notes)
for m = 1:N
  j = 0;
  lambda(m,1) = D(1,j+1)*exp(1i*2*pi*j*(m-1)/N);
  for j = 1:N-1
    lambda(m,1) = lambda(m,1) + D(1,j+1)*exp(1i*2*pi*j*(m-1)/N);
  end
  X(:,m) = exp(1i*(2*pi/L)*(m-1)*x);
end

% sanity check, X^-1*D*X should equal \Lambda
if norm(diag(lambda)-inv(X)*D*X) > 10*(N^2)*eps;
  error('Wrong eigenvalues/vectors!')
end

% initial condition
u_0 = real(X(:,2));
u_0 = exp(-((x-0.5)/0.1).^2);

% Solution in modal coordinates	
w_0 = inv(X)*u_0;
w = @(t,lambda,w_0) exp(-lambda.*t).*w_0;

% time interval 
t_max = 1;
N_t 	= 10;
dt 		= t_max/N_t;

figure
for n = 0:N_t	
	t = n*dt;
	u = X*w(t,lambda,w_0);	
	
	% plotting
	plot(x,u_0,'color',0.7*ones(3,1),'linewidth',4);
	hold on;
	
	% exact solution
  [xx,ii] = sort(mod(x+t,1));
  plot(xx,u_0(ii),'--','color',0*ones(3,1),'linewidth',2)

	plot(x,real(u),'color',[1 0 0],'linewidth',2);
	hold off;
	
	legend('IC','Exact sol.',sprintf('FD sol., t=%1.2f',t));
	xlabel('x'); ylabel('u');
    pause(.001);
	drawnow;
end 
