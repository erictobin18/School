function Final()
% I wrote a single finite difference solver and a single finite volume 
% solver that works for an n-dimensional problem by first order (global) 
% operator splitting

% Setup
clearvars; clc;
default_plot_attributes;
addpath('Euler_Equations');
addpath('Problem_Initialization');


nx = [50,100,150];
epsilon = [1,.5,.25];
gamma = 1.4;
fnum = 1; % figure number

problem1 = true;
problem2 = true;
problem3 = true;
problem4 = false;
 
%% Problem 1
if problem1 == true
for n = nx
    for e = epsilon
        mat(:,:,1) = -eye(3);
        mat(:,:,2) = eye(3);      
        bc.open_left  = {1, [0 1], mat};
        bc.open_right  = {n, [0 -1], mat};
        
        sol = Euler1D_FD(1,.2,n,@(x) Sod_Shock_Tube(x,gamma),bc,e,gamma);
        
        fnum = plot_1D(sol,fnum,false);
    end % for e
end % for n
end % if

%% Problem 2
if problem2 == true
for n = nx
    mat(:,:,1) = zeros(3); % dq = 0 (open)
    mat(:,:,2) = zeros(3);      
    bc.open_left  = {1, [0 1], mat};
    bc.open_right  = {n, [0 -1], mat};
    
    sol = Euler1D_FV(1,.2,n,@(x) Sod_Shock_Tube(x,gamma),bc,gamma);

    fnum = plot_1D(sol,fnum,false);
end % for n
end % if

%% Problem 3
if problem3 == true
for e = epsilon
    clear('mat');
    clear('bc');
    nx = 400;
    ny = 400;
    mat(:,:,1) = -eye(4);
    mat(:,:,2) = eye(4);
    bc.open_south = {[1:nx;ones(1,nx)], [0 0;0 1], mat};
    bc.open_west  = {[ones(1,ny);1:ny], [0 1;0 0], mat};
    bc.open_north = {[1:nx;ones(1,nx)*ny], [0 0;0 -1], mat};
    bc.open_east  = {[ones(1,ny)*nx;1:ny], [0 -1;0 0], mat};
    
    sol4_1 = Euler2D_FD(1,1,.3,nx,ny,@(x,y) Liska_Wendroff_4_1(x,y,gamma),bc,e,gamma);
%     sol4_2 = Euler2D_FD(1,1,.25,nx,ny,@(x,y) Liska_Wendroff_4_2(x,y,gamma),bc,e,gamma);
%     sol4_3 = Euler2D_FD(1,1,.2,nx,ny,@(x,y) Liska_Wendroff_4_3(x,y,gamma),bc,e,gamma);
%     sol_sr = Euler2D_FD(1,1,.2,nx,ny,@(x,y) Shock_Reflection(x,y,gamma),bc,e,gamma);

    fnum = plot_2D(sol4_1,fnum,true);
end % for e
end % if

%% Problem 4
if problem4 == true
nx = 400;
ny = 400;
mat(:,:,1) = zeros(4);
mat(:,:,2) = zeros(4);
bc.open_south = {[1:nx;ones(1,nx)], [0 0;0 1], mat};
bc.open_west  = {[ones(1,ny);1:ny], [0 1;0 0], mat};
bc.open_north = {[1:nx;ones(1,nx)*ny], [0 0;0 -1], mat};
bc.open_east  = {[ones(1,ny)*nx;1:ny], [0 -1;0 0], mat};

sol = Euler2D_FV(1,1,.2,nx,ny,@(x,y) Liska_Wendroff_4_1(x,y,gamma),bc);
% sol4_1 = Euler2D_FV(1,1,.3,nx,ny,@(x,y) Liska_Wendroff_4_1(x,y,gamma),bc);
% sol4_2 = Euler2D_FV(1,1,.25,nx,ny,@(x,y) Liska_Wendroff_4_2(x,y,gamma),bc);
% sol4_3 = Euler2D_FV(1,1,.2,nx,ny,@(x,y) Liska_Wendroff_4_3(x,y,gamma),bc);
% sol_sr = Euler2D_FD(1,1,.2,nx,ny,@(x,y) Shock_Reflection(x,y,gamma),bc,e,gamma);

fnum = plot_2D(sol,fnum,true);
end % if

end % final

%% Wrapper Functions
% 1(a)
function solution = Euler1D_FD(x, t, nx, q0fun, bc, epsilon,gamma) % wrapper function
temp = linspace(0,x,nx + 1);
xvect = temp(1:end-1);
f = @(q,t) EulerFun1D(q,t,gamma);
[tsol, qsol] = solve_FD({xvect'},[0 t/100 t],{f},q0fun(xvect),epsilon,bc);
solution = {tsol,xvect,qsol};
end
% 3(a)
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
[tsol, qsol] = solve_FV({xvect'},[0 t/100 t],{f},q0fun(xvect),bc);
solution = {tsol,xvect,qsol};
end
% 4(a)
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

%% Plotting
function [fnum] = plot_1D(sol,fnum,animated)
%% exact
% exact.txt downloaded and stripped of its first row to eliminate column
% headings
ex = importdata('exact.txt');
gamma = 1.4;

exx = ex(:,1);
exrho = ex(:,2);
exp = ex(:,3);
exu = ex(:,4);
exe = exp./(exrho*(gamma - 1));

ts = sol{1,1};
xs = sol{1,2};
qs = sol{1,3};

%% Plotting
figure(fnum);
fnum = fnum + 1;

start = 1;
if ~animated
    start = length(ts);
end
for i=start:length(ts)
    rho = permute(qs(1,:,i),[2,1,3]); % Conservative variables
    momx = permute(qs(2,:,i),[2,1,3]);
    rho_eng = permute(qs(3,:,i),[2,1,3]);

    u = momx./rho;
    eng = rho_eng./rho;
  subplot(1,3,1); 
    plot(xs,rho,'k'); ylim([0 1]); xlabel('x'); ylabel('rho');
    hold on;
    plot(exx,exrho,'--k');
    hold off;
  subplot(1,3,2); 
    plot(xs,u,'k'); ylim([0 1.2]); xlabel('x'); ylabel('u');
    hold on;
    plot(exx,exu,'--k');
    hold off;
  subplot(1,3,3);
    plot(xs,eng,'k'); ylim([0 5]); xlabel('x'); ylabel('energy');
    hold on;
    plot(exx,exe,'--k');
    hold off;
  pause(0.001);
    drawnow;
end % for i
end

function [fnum] = plot_2D(sol,fnum,animated)
ts = sol{1,1};
xs = sol{1,2};
ys = sol{1,3};
qs = sol{1,4};

%% Plotting
figure(fnum);
fnum = fnum + 1;

start = 1;
if ~animated
    start = length(ts);
end
for i=start:length(ts)
    rho = permute(qs(1,:,:,i),[2,3,1,4]); % Conservative variables
    momx = permute(qs(2,:,:,i),[2,3,1,4]);
    momy = permute(qs(3,:,:,i),[2,3,1,4]);
    eng = permute(qs(4,:,:,i),[2,3,1,4]);

    subplot(2,2,3); 
    mesh(xs,ys,rho); zlim([0 1]); xlabel('x'); ylabel('y'); zlabel('rho');
%         hold on;
%         plot(exx,exrho,'--k');
%         hold off;
  subplot(2,2,1); 
    mesh(xs,ys,momx); zlim([0 3]); xlabel('x'); ylabel('y'); zlabel('mom_x');
%         hold on;
%         plot(exx,exu,'--k');
%         hold off;
  subplot(2,2,2);
    mesh(xs,ys,momy); zlim([0 3]); xlabel('x'); ylabel('y'); zlabel('mom_y');
%         hold on;
%         plot(exx,exe,'--k');
%         hold off;
  subplot(2,2,4);
    mesh(xs,ys,eng); zlim([0 5]); xlabel('x'); ylabel('y'); zlabel('energy');
%         hold on;
%         plot(exx,exe,'--k');
%         hold off;

  pause(0.001);
    drawnow;
end % for i
end
