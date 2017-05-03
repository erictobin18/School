clear all; clc;

%=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Burger's equation
%
%   du/dt = -0.5*d/dx(u^2) , x \in [0,1)
%
% with periodic boundary conditions. 
%
%=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------
  % spatial discretization
	N   = 201;
	L   = 1;
	dx  = L/N;
	x   = dx*linspace(0,N-1,N).';

  % time interval 
  t_max  = 1;
  N_t    = 200;
  dt     = t_max/N_t;
  t_plot = [0:dt:t_max];

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
D = D_upwind_3rd;

% initial condition
u_0 = 1 + exp(-((x-0.5)/0.1).^2);

% time integration   
[t,u_con]     = ode45(@(t,u) lecture7_con_burgers_ODE(t,u,D),t_plot,u_0); 
[t,u_noncon]  = ode45(@(t,u) lecture7_noncon_burgers_ODE(t,u,D),t_plot,u_0); 


figure;
for n = 1:length(t)     
  
  % Burger's equation solution
  plot(x,u_con(n,:),    'color',[1 0 0],'linewidth',2);
  hold on;
  plot(x,u_noncon(n,:),  'color',[0 0 1],'linewidth',2);
  hold off;
  
  legend( 'Con. Burgers equation solution',...
      'Non-Con. Burgers equation solution',...
        'location','southoutside','orientation','vertical');

  xlabel('x'); ylabel('u');
  xlim([0 1]); ylim([0.9 2.1]);
  
  drawnow;
end
