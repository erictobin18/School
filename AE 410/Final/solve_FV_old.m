function [tsol, ysol] = solve_FV_old(x, t, flux, y0, bc) % implementation space
dt = .002;
tsol = t(1):dt:t(2); % Guess what the timesteps will be
tc = t(1);

% x = x(:); % ensure x is col vector
% dx = circshift(x,-1) - x; % vector of delta x for interval to the right of the points
% dx(end) = mean(dx(1:end-1)); % imaginary extra delta x

n = length(x);

yc = y0;
ycg = y0;
ysol = zeros(size(yc,1),size(yc,2),length(tsol)); % at each timestep, save yc as a row vector

iter = 1; % iteration number
while tc < t(2) + dt
    ysol(:,:,iter) = yc; % save yc
    tsol(iter) = tc; % save time of current iteration
    
    f = zeros(3,n);
    fg = zeros(3,n);
    if strcmp(bc, 'periodic')
%         f = flux(yc,circshift(yc,[0,-1]),t);
        for i=1:n-1
            fg(:,i) = AppRiemannSol(ycg(:,i),ycg(:,i+1));
        end
        fg(:,n) = AppRiemannSol(ycg(:,n),ycg(:,1));
    end
    
%     yn = yc + dt*(circshift(f,[0,1]) - f);
    yng = ycg + dt*(circshift(fg,[0,1]) - fg);
    % update "current" state yc and time tc
    tc = tc + dt;
%     yc = yn;
    ycg = yng;
    yc = ycg;
    iter = iter + 1;
    if tc >= t(2) - dt
        disp('hello');
    end
end

end

function ynew = characteristic(amat, yc, d, dt)
eigvect = zeros(3,3,n);
eigval = zeros(3,3,n);
for i=1:n
    [eigvect(:,:,i),eigval(:,:,i)] = eig(amat(i));
end

% Characteristic variables
wc = zeros(3,n);
for i=1:n
    wc(:,i) = eigvect(:,:,i)\yc(:,i);
end
wn = zeros(3,n);
for i=1:3
    lambda = zeros(n);
    for j=1:n
        lambda(j,j) = eigval(i,i,j);
    end
    wn(i,:) = expm(-lambda*d*dt)*wc(i,:);
end

% Primitive variables
ynew = zeros(3,n);
for i=1:n
    ynew(:,i) = eigvect(:,:,i)*wn(:,i);
end
end
