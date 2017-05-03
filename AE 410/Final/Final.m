function Final()
% I wrote a single finite difference solver that works for an n-dimensional
% problem by first order (global) operator splitting

% Setup

clearvars; clc;
default_plot_attributes;
addpath('Euler_Equations');


nx = [50,100,150];
epsilon = [1,.5,.25];
gamma = 1.4;
fnum = 1; % figure number

problem1 = false;
problem2 = false;
problem3 = false;
problem4 = true;


%% exact
% exact.txt downloaded and stripped of its first row to eliminate column
% headings
ex = importdata('exact.txt');

exx = ex(:,1);
exrho = ex(:,2);
exp = ex(:,3);
exu = ex(:,4);
exe = exp./(exrho*(gamma - 1));
 
%% Problem 1
if problem1 == true
for n = nx
    for e = epsilon
        mat(:,:,1) = -eye(3);
        mat(:,:,2) = eye(3);      
        bc.open_left  = {[1], [0 1], mat};
        bc.open_right  = {[n], [0 -1], mat};
        
        sol = Euler1D_FD(1,.2,n,@(x) Euler1DInit(x,gamma),bc,e,gamma);
        
        ts = sol{1,1};
        xs = sol{1,2};
        qs = sol{1,3};

        % Plotting
        figure(fnum);
        fnum = fnum + 1;
        for i = 1:length(ts)	
          u = qs(1,:,i); % rho
          v = qs(2,:,i); % rho*u
          w = qs(3,:,i); % rho*E


          subplot(1,3,1); 
            plot(xs,u,'k'); ylim([0 1]); xlabel('x'); ylabel('rho');
            hold on;
            plot(exx,exrho,'--k');
            hold off;
          subplot(1,3,2); 
            plot(xs,v./u); ylim([-1 1]); xlabel('x'); ylabel('u');
            hold on;
            plot(exx,exu,'--k');
            hold off;
          subplot(1,3,3);
            plot(xs,w./u); ylim([0 3]); xlabel('x'); ylabel('E');
            hold on;
            plot(exx,exe,'--k');
            hold off;

          pause(0.001);

            drawnow;
        end % for i
    end % for e
end % for n
end % if

%% Problem 2
if problem2 == true
for n = nx
    mat(:,:,1) = zeros(3); % dq = 0 (open)
    mat(:,:,2) = zeros(3);      
    bc.open_left  = {[1], [0 1], mat};
    bc.open_right  = {[n], [0 -1], mat};
    
    sol = Euler1D_FV(1,.2,n,@(x) Euler1DInit(x,gamma),bc,gamma);

    ts = sol{1,1};
    xs = sol{1,2};
    qs = sol{1,3};

    %% Plotting
    figure(fnum);
    fnum = fnum + 1;
    for i = 1:length(ts)	
      u = qs(1,:,i); % rho
      v = qs(2,:,i); % rho*u
      w = qs(3,:,i); % rho*E


      subplot(1,3,1); 
        plot(xs,u,'k'); ylim([0 1]); xlabel('x'); ylabel('rho');
        hold on;
        plot(exx,exrho,'--k');
        hold off;
      subplot(1,3,2); 
        plot(xs,v./u); ylim([-1 1]); xlabel('x'); ylabel('u');
        hold on;
        plot(exx,exu,'--k');
        hold off;
      subplot(1,3,3);
        plot(xs,w./u); ylim([0 3]); xlabel('x'); ylabel('E');
        hold on;
        plot(exx,exe,'--k');
        hold off;

      pause(0.001);

        drawnow;
    end % for i
end % for n
end % if

%% Problem 3
if problem3 == true
for e = epsilon
    nx = 400;
    ny = 400;
    mat(:,:,1) = -eye(4);
    mat(:,:,2) = eye(4);
    bc.open_south = {[1:nx;ones(1,nx)], [0 0;0 1], mat};
    bc.open_west  = {[ones(1,ny);1:ny], [0 1;0 0], mat};
    bc.open_north = {[1:nx;ones(1,nx)*ny], [0 0;0 -1], mat};
    bc.open_east  = {[ones(1,ny)*nx;1:ny], [0 -1;0 0], mat};
    
    sol = Euler2D_FD(1,1,.2,nx,ny,@(x,y) Euler2DInit(x,y,gamma),bc,e,gamma);
    
    ts = sol{1,1};
    xs = sol{1,2};
    qs = sol{1,3};
    qs = sol{1,4};
    
    %% Plotting
    
    figure(fnum);
    fnum = fnum + 1;
    for i=1:length(ts)
        rho = qs(1,:,:,i); % Conservative variables
        momx = qs(2,:,:,i);
        momy = qs(3,:,:,i);
        eng = qs(4,:,:,i);
        
        
        subplot(2,2,3); 
        mesh(xs,qs,rho); xlabel('x'); ylabel('y'); zlabel('rho');
%         hold on;
%         plot(exx,exrho,'--k');
%         hold off;
      subplot(2,2,1); 
        plot(xs,qs,momx); xlabel('x'); ylabel('y'); zlabel('mom_x');
%         hold on;
%         plot(exx,exu,'--k');
%         hold off;
      subplot(2,2,2);
        plot(xs,qs,momy); xlabel('x'); ylabel('y'); zlabel('mom_y');
%         hold on;
%         plot(exx,exe,'--k');
%         hold off;
      subplot(2,2,4);
        plot(xs,qs,eng); ylim([0 3]); xlabel('x'); ylabel('y'); zlabel('energy');
%         hold on;
%         plot(exx,exe,'--k');
%         hold off;

      pause(0.001);

        drawnow;
    end % for i
end % for e
end % if

%% Problem 4
if problem4 == true
nx = 400;
ny = 400;
mat(:,:,1) = -eye(4);
mat(:,:,2) = eye(4);
bc.open_south = {[1:nx;ones(1,nx)], [0 0;0 1], mat};
bc.open_west  = {[ones(1,ny);1:ny], [0 1;0 0], mat};
bc.open_north = {[1:nx;ones(1,nx)*ny], [0 0;0 -1], mat};
bc.open_east  = {[ones(1,ny)*nx;1:ny], [0 -1;0 0], mat};

sol = Euler2D_FV(1,1,.2,nx,ny,@(x,y) Euler2DInit(x,y,gamma),bc);

ts = sol{1,1};
xs = sol{1,2};
qs = sol{1,3};
qs = sol{1,4};

%% Plotting

figure(fnum);
fnum = fnum + 1;
for i=1:length(ts)
    rho = qs(1,:,:,i); % Conservative variables
    momx = qs(2,:,:,i);
    momy = qs(3,:,:,i);
    eng = qs(4,:,:,i);


    subplot(2,2,3); 
    mesh(xs,qs,rho); xlabel('x'); ylabel('y'); zlabel('rho');
%         hold on;
%         plot(exx,exrho,'--k');
%         hold off;
  subplot(2,2,1); 
    plot(xs,qs,momx); xlabel('x'); ylabel('y'); zlabel('mom_x');
%         hold on;
%         plot(exx,exu,'--k');
%         hold off;
  subplot(2,2,2);
    plot(xs,qs,momy); xlabel('x'); ylabel('y'); zlabel('mom_y');
%         hold on;
%         plot(exx,exe,'--k');
%         hold off;
  subplot(2,2,4);
    plot(xs,qs,eng); ylim([0 3]); xlabel('x'); ylabel('y'); zlabel('energy');
%         hold on;
%         plot(exx,exe,'--k');
%         hold off;

  pause(0.001);

    drawnow;
end % for i
end % if

end % final
    
% 1(a)
function solution = Euler1D_FD(x, t, nx, q0fun, bc, epsilon,gamma) % wrapper function
temp = linspace(0,x,nx + 1);
xvect = temp(1:end-1);
f = @(q,t) EulerFun1D(q,t,gamma);
[tsol, qsol] = solve_FD({xvect'},[0 t/100 t],{f},q0fun(xvect),epsilon,bc);
solution = {tsol,xvect,qsol};
end
% 3
function solution = Euler2D_FD(x, y, t, nx, ny, q0fun, bc, epsilon, gamma) % wrapper function
tempx = linspace(0,x,nx + 1);
xvect = tempx(1:end-1);
tempy = linspace(0,y,ny + 1);
yvect = tempy(1:end-1);
f = @(q,t) EulerFunX2D(q,t,gamma);
g = @(q,t) EulerFunY2D(q,t,gamma);

[tsol, qsol] = solve_FD({xvect',yvect'},[0 t/100 t],{f,g},q0fun(xvect,yvect),epsilon,bc);
solution = {tsol,xvect,yvect,qsol};
end
% 2
function solution = Euler1D_FV(x, t, nx, q0fun, bc, gamma) % wrapper function
temp = linspace(0,x,nx + 1);
xvect = temp(1:end-1);
f = @(q_l,q_r,t) Flux_Roe(q_l,q_r,t,gamma);
[tsol, qsol] = solve_FV({xvect'},[0 t/169 t],{f},q0fun(xvect),bc);
solution = {tsol,xvect,qsol};
end
% 4
function solution = Euler2D_FV(x, y, t, nx, ny, q0fun, bc) % wrapper function
tempx = linspace(0,x,nx + 1);
xvect = tempx(1:end-1);
tempy = linspace(0,y,ny + 1);
yvect = tempy(1:end-1);
f = @(q_l,q_r,t) RoeFlux_for_Q4_wrapperX2D(q_l,q_r);
g = @(q_l,q_r,t) RoeFlux_for_Q4_wrapperY2D(q_l,q_r);

[tsol, qsol] = solve_FV({xvect',yvect'},[0 t/100 t],{f,g},q0fun(xvect,yvect),bc);
solution = {tsol,xvect,yvect,qsol};
end

function flux = RoeFlux_for_Q4_wrapperX2D(q_l,q_r)
flux = zeros(size(q_l));
q_l(:,:,5) = q_l(:,:,4);
q_l(:,:,4) = zeros(size(q_l(:,:,4)));
q_r(:,:,5) = q_r(:,:,4);
q_r(:,:,4) = zeros(size(q_r(:,:,4)));

for i=1:size(q_l,1)
    for j=1:size(q_l,2)
        left = q_l(i,j,:);
        right = q_r(i,j,:);
        temp = RoeFlux_for_Q4(left(:),right(:))';
        flux(i,j,:) = [temp(1:3),temp(5)];
    end
end
end

function flux = RoeFlux_for_Q4_wrapperY2D(q_l,q_r) % flux in y direction!
% swap x and y, load into RoeFlux, then swap x and y again
flux = zeros(size(q_l));
q_l(:,:,5) = q_l(:,:,4); % 5 <= 4
q_l(:,:,4) = q_l(:,:,2); % 4 <= 2
q_l(:,:,2) = q_l(:,:,3); % 2 <= 3
q_l(:,:,3) = q_l(:,:,4); % 3 <= 4
zeros(size(q_l(:,:,4))); % 4 <= zeros
q_r(:,:,5) = q_r(:,:,4); % 5 <= 4
q_r(:,:,4) = q_r(:,:,2); % 4 <= 2
q_r(:,:,2) = q_r(:,:,3); % 2 <= 3
q_r(:,:,3) = q_r(:,:,4); % 3 <= 4
zeros(size(q_r(:,:,4))); % 4 <= zeros

for i=1:size(q_l,1)
    for j=1:size(q_l,2)
        left = q_l(i,j,:);
        right = q_r(i,j,:);
        temp = RoeFlux_for_Q4(left(:),right(:))';
        flux(i,j,:) = [temp(1),temp(3),temp(2),temp(5)];
    end
end
end

function out = Euler1DInit(x,gamma) % initial conditions function
% variables starting with "y" always refers to a Nx3 matrix
% the three columns are columns for rho, rho*u, rho*E
% in other words, instead of something like 
% y(1:N) = rho; y(N+1,2*N) = rho_u; y(2*N+1,3*N) = rho_E;
% this program arranges the data like
% y = [rho, rho_u, rho_E]
% as a matrix. On disk, they are stored the same; the only difference is
% how it is indexed
n = length(x);
n2 = floor(n/2);
rho1 = 1;
p1 = 1;
u1 = 0;
rho2 = .125;
p2 = .1;
u2 = 0;
out = [ 
    ones(n2,1)*rho1   ones(n2,1)*rho1*u1   ones(n2,1)*ideal_rho_E(p1,rho1,rho1*u1,gamma);
    ones(n-n2,1)*rho2 ones(n-n2,1)*rho2*u2 ones(n-n2,1)*ideal_rho_E(p2,rho2,rho2*u2,gamma);
    ]';
end

function out = Euler2DInit(x,y,gamma) % initial conditions function
grid = ones(4,length(x),length(y));
out = grid;
end

function out = pack(t,y,gamma,nx)

temp = Flux_Roe(...
    [...
        y(1:nx)';
        y(nx+1:2*nx)';
        y(2*nx+1:3*nx)';
        ],...
    circshift([
        y(1:nx)';
        y(nx+1:2*nx)';
        y(2*nx+1:3*nx)';
        ],[0,-1]),...
    t,...
    gamma);

out = [temp(1,:)';temp(2,:)';temp(3,:)'];
    
end

function out = enthalpy(y,g)
out = y(3,:)./y(1,:) - (1/2)*(y(2,:)./y(1,:)).^2 + ideal_p(y,g)./y(1,:);
end
function out = ideal_p(y,g) % solve equation of state for p
out = (g - 1)*(y(3,:) - (1/2)*(y(2,:).^2)./y(1,:));
end
function out = ideal_rho_E(p,rho,rho_u,gamma) % solve equation of state for rho*E
out = (p/(gamma - 1)) + (1/2)*(rho_u.^2)./rho;
end
