clear all;
clc;

% -------------------------------------------
  addpath ./matrix_operators
  
  % set default plot attributes
   default_plot_attributes

	Nx = 100;
	L  = 1;
	dx = L/Nx;
	x  = dx*linspace(0,Nx-1,Nx).';
	
% --------------------------------------------

% %%%%%%%%%% Question 2 %%%%%%%%%%%%%%%%%%%%%%

% construct periodic FD matrix
D{1} = upwind_1st_order(x,1);
D{2} = upwind_2nd_order(x,1);
D{3} = central_2nd_order(x,1);
D{4} = upwind_3rd_order(x,1);

% initial condition
u_0(x>1/3,1) = 1;
u_0(x>2/3,1) = 0;

% exact solution
figure; 
plot(x,u_0,'color',0.7*ones(3,1),'linewidth',2);
hold on;

for i = 1:4
	[X,Lambda] = eig(full(-D{i}));
	lambda = diag(Lambda);
		
	w_0 = inv(X)*u_0;
	w = @(t,lambda,w_0) exp(lambda.*t).*w_0;
	u = X*w(1,lambda,w_0);	
	plot(x,real(u),'linewidth',1);
end

xlabel('x'); ylabel('u');
legend(	'exact',...
		'1^{st}-order upwind',...
		'2^{nd}-order upwind',...
		'2^{nd}-order central',...
		'3^{rd}-order upwind'); 

% %%%%%%%%%% Question 3 %%%%%%%%%%%%%%%%%%%%%%
y = x;
for i = 1:Nx
	for j = 1:Nx
		U(i,j) = exp(-((x(i)-0.45).^2/0.05^2+(y(j)-0.4).^2/0.15^2));
	end
end

D = D{3};

dfdx 	= D*U;
dfdy 	= U*D';
d2fdxy  = D*U*D';

figure;
subplot(2,2,1); imagesc(x,y,U); 
axis square; shading interp; title('f(x,y)'); xlabel('y'); ylabel('x');
subplot(2,2,2); imagesc(x,y,dfdx); 
axis square; shading interp; title('df/dx'); xlabel('y'); ylabel('x'); 
subplot(2,2,3); imagesc(x,y,dfdy); 
axis square; shading interp; title('df/dy'); xlabel('y'); ylabel('x'); 
subplot(2,2,4); imagesc(x,y,d2fdxy); 
axis square; shading interp; title('d^2f/dyx'); xlabel('y'); ylabel('x'); 
