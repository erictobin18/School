function HW4
clearvars; clc;
global gamma
gamma = 1.4;
% Copied from HW Code:
% -------------------------------------------
  % include path for matrix operator m files
  addpath ./matrix_operators

  % set default plot attributes
  default_plot_attributes  

% exact.txt downloaded and stripped of its first row to eliminate column
% headings
ex = importdata('exact.txt');

exx = ex(:,1);
exrho = ex(:,2);
exp = ex(:,3);
exu = ex(:,4);
exe = exp./(exrho*(gamma - 1));

x = (0:.005:1);
%x = (0:.1:1); % run coarsely discretized system for simpler debugging, not 
               % stable even for Riemann Solvers
deltax = x(2)-x(1);
n = length(x);

% Initial Conditions
rho1 = 1;
p1 = 1;
u1 = 0;
mom1 = rho1*u1;
eng1 = p1/(gamma-1) + (1/2)*rho1*u1^2;

rho2 = .125;
p2 = .1;
u2 = 0;
mom2 = rho2*u2;
eng2 = p2/(gamma-1) + (1/2)*rho2*u2^2;

rho = [ones((n+1)/2,1)*rho1;ones((n-1)/2,1)*rho2];
mom = [ones((n+1)/2,1)*mom1;ones((n-1)/2,1)*mom2];
eng = [ones((n+1)/2,1)*eng1;ones((n-1)/2,1)*eng2];

v = [rho; mom; eng]; % initial condition

% 1(a)
[~,y1] = ode45(@(t,v) get_derivative(t,v,deltax,n,@avg), [0 0.2], v);
% 1(b)
[~,y2] = ode45(@(t,v) get_derivative(t,v,deltax,n,@riemann_fixed), [0 0.002], v);
% 2
[~,y3] = ode45(@(t,v) get_derivative(t,v,deltax,n,@riemann_wall), [0 0.5], v);

% unpack data
rho1 = y1(end,1:n);
rho2 = y2(end,1:n);
rho3 = y3(end,1:n);

mom1 = y1(end,n+1:2*n);
mom2 = y2(end,n+1:2*n);
mom3 = y3(end,n+1:2*n);

eng1 = y1(end,2*n+1:3*n);
eng2 = y2(end,2*n+1:3*n);
eng3 = y3(end,2*n+1:3*n);

% compute primitive (non-conservative) variables
u1 = mom1./rho1;
u2 = mom2./rho2;
u3 = mom3./rho3;

p1 = (gamma - 1)*(eng1 - (1/2)*mom1.*rho1);
p2 = (gamma - 1)*(eng2 - (1/2)*mom2.*rho2);
p3 = (gamma - 1)*(eng3 - (1/2)*mom3.*rho3);

e1 = p1./(rho1*(gamma - 1));
e2 = p2./(rho2*(gamma - 1));
e3 = p3./(rho3*(gamma - 1));

% AVERAGING
figure(1);
subplot(2,1,1);
plot(x,rho1,'k');
grid on;
xlabel('x');
ylabel('rho');
subplot(2,1,2);
plot(x,u1,'k');
grid on;
xlabel('x');
ylabel('u');
ylim([-1,1]);

figure(2);
subplot(2,1,1);
plot(x,p1,'k');
grid on;
xlabel('x');
ylabel('p');
ylim([0,1]);
subplot(2,1,2);
plot(x,e1,'k');
grid on;
xlabel('x');
ylabel('e');
ylim([0,4]);

% FIXED BOUNDARY CONDITIONS
figure(3);
subplot(2,1,1);
hold on;
plot(x,rho2,'k','LineWidth',1);
plot(exx,exrho,'k--','LineWidth',2);
hold off;
grid on;
xlabel('x');
ylabel('rho');
legend('Computed','Exact');
subplot(2,1,2);
hold on;
plot(x,u2,'k','LineWidth',1);
plot(exx,exu,'k--','LineWidth',2);
hold off;
grid on;
xlabel('x');
ylabel('u');
legend('Computed','Exact');

figure(4);
subplot(2,1,1);
hold on;
plot(x,p2,'k','LineWidth',1);
plot(exx,exp,'k--','LineWidth',2);
hold off;
grid on;
xlabel('x');
ylabel('p');
legend('Computed','Exact');
subplot(2,1,2);
hold on;
plot(x,e2,'k','LineWidth',1);
plot(exx,exe,'k--','LineWidth',2);
hold off;
grid on;
xlabel('x');
ylabel('e');
legend('Computed','Exact');

% WALL BOUNDARY CONDITIONS
figure(5);
subplot(2,1,1);
plot(x,rho3,'k');
grid on;
xlabel('x');
ylabel('rho');
subplot(2,1,2);
plot(x,u3,'k');
grid on;
xlabel('x');
ylabel('u');

figure(6);
subplot(2,1,1);
plot(x,p3,'k');
grid on;
xlabel('x');
ylabel('p');
subplot(2,1,2);
plot(x,e3,'k');
grid on;
xlabel('x');
ylabel('e');

% Uncomment for animated plots
% Copied from Homework Code
% for i = 1:length(t)	
% 
%   u_1 = y(i,1:n);
%   u_2 = y(i,n+1:2*n);
%   u_3 = y(i,2*n+1:3*n);
% 
%   subplot(1,3,1); 
%     plot(x,u_1); ylim([-1 1]); xlabel('x'); ylabel('rho');
%   subplot(1,3,2); 
%     plot(x,u_2); ylim([-1 1]); xlabel('x'); ylabel('momentum');
%   subplot(1,3,3); 
%     plot(x,u_3); ylim([-1 1]); xlabel('x'); ylabel('energy');
%   
%   pause(0.001);
% 
% 	drawnow;
% end 


end

function vdot = get_derivative(~,v,deltax,n,func)
% Unpack
rho = v(1:n);
mom = v(n + 1:2*n);
eng = v(2*n + 1:3*n);

net_flux = func(rho,mom,eng,n);

vdot = net_flux/deltax;
end

function out = avg(rho, mom, eng, n)
global gamma

% Compute F at interval centers (average)
f_rho = mom; % F_1 = q2
f_mom = (mom.^2)./rho + (gamma - 1).*(eng - (1/2)*(mom.^2)./rho); % F_2 = q2^2/q1 + (g-1)(q3-(1/2)q2^2/q1)
f_eng = ((mom./rho).*(eng + (gamma - 1)*(eng - (1/2)*(mom.^2)./rho))); % F_3 = (q2/q1)(q3 + (g-1)(q3-(1/2)q2^2/q1))

% Use F in adjacent cells (similar to central second order scheme)
fp = [circshift(f_rho,+1); circshift(f_mom,+1); circshift(f_eng,+1)];
fm = [circshift(f_rho,-1); circshift(f_mom,-1); circshift(f_eng,-1)];

out = (fp - fm)/2; % Identical to second order FD

out(1) = 0; out(n) = 0; out(n+1) = 0; out(2*n) = 0; out(2*n+1) = 0; out(3*n)=0; % Fixed BC
end

function out = riemann_fixed(rho, mom, eng, n)
% Find the interface (F) values (periodic boundary)
f_rho = zeros(n - 1,1);
f_mom = zeros(n - 1,1);
f_eng = zeros(n - 1,1);

for i = 1:n-1 % there are N-1 interior F values
    temp = AppRiemannSol(...
        [rho(i);mom(i);eng(i)],...
        [rho(i+1);mom(i+1);eng(i+1)]);
    f_rho(i) = temp(1);
    f_mom(i) = temp(2);
    f_eng(i) = temp(3);
end

% Offset values to the left or right to get F_L and F_R
fp = [f_rho; 0; f_mom; 0; f_eng; 0];
fm = [0; f_rho; 0; f_mom; 0; f_eng];

out = fm - fp;

out(1) = 0; out(n) = 0; out(n+1) = 0; out(2*n) = 0; out(2*n+1) = 0; out(3*n)=0; % Fixed BC
end

function out = riemann_wall(rho, mom, eng, n)
% Identical to riemann_fixed except last lines (BC)
% Find the interface (F) values (periodic boundary)
f_rho = zeros(n - 1,1);
f_mom = zeros(n - 1,1);
f_eng = zeros(n - 1,1);

for i = 1:n-1 % there are N-1 interior F values
    temp = AppRiemannSol(...
        [rho(i);mom(i);eng(i)],...
        [rho(i+1);mom(i+1);eng(i+1)]);
    f_rho(i) = temp(1);
    f_mom(i) = temp(2);
    f_eng(i) = temp(3);
end

fp = [f_rho; 0; f_mom; 0; f_eng; 0];
fm = [0; f_rho; 0; f_mom; 0; f_eng];

out = fm - fp;

% Wall Boundary Conditions
out(1) = out(2); out(n) = out(n-1); % Reflect Density
out(n+1) = -out(n+2); out(2*n) = -out(2*n-1); % Negative Reflect Momentum
out(2*n+1) = out(2*n+2); out(3*n)= out(3*n-1); % Reflect Energy
end