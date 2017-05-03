function dqdt = problem3_system(~,q,N,gamma,Dp,Dm)

rho = q(1  :  N);
u = q(N+1  :2*N); 
p = q(2*N+1:3*N);

drp = Dp*rho;
drm = Dm*rho;
dup = Dp*u;
dum = Dm*u;
dpp = Dp*p;
dpm = Dm*p;

c = sqrt(gamma*p./rho);

dqdt = zeros(3*N,1);

% Lt = zeros(N,N,N);
% Lt(:,1,1) = u;
% Lt(:,2,2) = u-c;
% Lt(:,3,3) = u+c;
% 
% Lp = Lt.*logical(Lt>=0);
% Lm = Lt.*logical(Lt<0);
% 

%X(m) = [1,1./c(m).^2,1./c(m).^2;0,1./(c(m).*rho(m)),-1./(c(m).*rho(m));0,1,1];

%dwp(m) = X(m)\[drp(m);dup(m);dpp(m)];

for m=1:N
    L = diag([u(m),u(m)-c(m),u(m)+c(m)]);
    X = [1,1/c(m)^2,1/c(m)^2;0,1/(c(m)*rho(m)),-1/(c(m)*rho(m));0,1,1];
    Lp = L.*logical(L>=0);
    Lm = L.*logical(L<0);
    dwp = X\[drp(m);dup(m);dpp(m)];
    dwm = X\[drm(m);dum(m);dpm(m)];
    dw = Lp*dwp + Lm*dwm;
    dq = X*dw;
    dqdt(m) = dq(1);
    dqdt(m+N) = dq(2);
    dqdt(m+2*N) = dq(3); 
end

end


