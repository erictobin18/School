function [tsol,qsol] = solve_FD(grid, t, f, q0, epsilon, bc) % implementation space
%SOLVE_FD Solves a finite difference problem.
% grid is expected to be a cell array of column vectors representing the
% gridpoints of the simulation. Grids must be rectangular, but can be
% unevenly spaced.

% t is expected to be a matrix with two or three values. The first value is
% the start time of the simulation, and the last is the end time. If a
% third value is present between these, it is the timestep dt.

% f is expected to be either a cell array of function handles or a single
% function handle. If it is an array, each handle corresponds to one
% dimension of the grid. If it is a single handle, the single handle is
% applied to all dimensions of the grid. The function args are expected to
% be the current state, q, and the current time. The current state is
% stored with the the solution variables along the first dimension of a
% matrix. This state is evaluated over the entire grid, with each dimension
% of the grid corresponding to an additional dimension of q. So, a 2D grid
% of 400 points with 4 solution variables stores the currents state in a
% 4x400x400 matrix.

% q0 is expected to be the current state at time = t(1). Thus, it must be
% the same size as q.

% epsilon is expected to be either a single number or a vector containing
% the artificial viscosity coefficients for each dimension

% bc is expected to be a structure. Each field is a boundary condition,
% whose value is a cell array. The first value of the array is a matrix of
% coordinates of the nodes that the bc is applied to. The second value is a
% matrix indicating the "stencil" for finding the time derivative of the
% state of each node in the bc. Each column in this matrix gives the
% coordinates relative to the bc node of a neighbor node. The last value is
% a vector of matrices indicating the weight each neighboring node should
% have in computing the derivative. The vector corresponds to the list of
% neighbors, and the weight is found by multiplying the matrix for each
% neighbor by the state at the neighbor.
% Check inputs

%error('add type checking with inputparser');
disp('WARN: solve_FD not checking types');



% Setup
dim = length(grid);
g_shape = repmat({':'},1,dim); % g_shape{:} will index a matrix to return all dim dimensions
n_eq = size(q0,1); % Number of variables to solve for

n = zeros(dim,1);
for i=1:dim
    n(i) = length(grid{i});
end

dt = .001;
tc = t(1);
tf = t(end);
if length(t) == 3
    dt = t(2);
elseif length(t) == 1
    tc = 0;
elseif length(t) ~= 2
    error('Incorrect Length "t" vector');
end
tsol = tc:dt:tf; % Guess what the timesteps will be
nt = length(tsol);

dgr = grid; % space between adjacent grid points to the "right" (i,i+1)
dgl = grid; % space between adjacent grid points to the "left" (i-1,i)
dga = grid; % average of dgr and dgl = (i+1 + i-1)/2
for i=1:dim
    temp = grid{i};
    dgr{i} = circshift(temp,-1) - temp; % vector of delta x for interval to the right of the points
    dgl{i} = temp - circshift(temp,1);
    dgr{i}(end) = mean(dgr{i}(1:end-1)); % imaginary extra delta x
    dgl{i}(1) = mean(dgl{i}(2:end));
    dga{i} = (dgr{i} + dgl{i})/2;
end

qc = q0;
qsol = zeros([n_eq, n', nt]); % save each variable at every grid point at every timestep

d = cell(1,dim);
d2 = cell(1,dim);
for i=1:dim
    d{i} = central_fd(grid{i},'periodic');
    d2{i} = central_fd_second(grid{i},'periodic');
end

iter = 1; % iteration number
while tc < tf + dt
    if mod(iter,10)==0
        disp(strcat('Iteration ',num2str(iter)));
    end
    qsol(:,g_shape{:},iter) = qc; % save qc
    tsol(iter) = tc; % save time of current iteration
    
    for i=1:dim % Lie Splitting: Solve each dimension independently first-order
        % transpose the current values so the dimension of interest is
        % along the rows (dimension 1 is for the solution vars):
        rot = circshift((1:dim)',1-i)';
        unrot = circshift((1:dim)',i-1)';
        qtemp = permute(qc,[1, rot+1]);
                         % col % row % the rest
        fvects = f{i}(qtemp,t); % array of F(q) vectors at each point
        
        % Transpose fvect so rows (dim 2) are cols (dim 1) and find spatial derivative
        temp = permute(fvects,[2:dim+1,1]);
        s = size(temp);
        dfdx = d{i}*reshape(temp,n(i),[]);
                                   % circular shift of dim+1 indices
        clear('temp');
        % Transpose q so rows are cols and find second spatial derivatives
        d2qdx2 = d2{i}*reshape(permute(qtemp, [2:dim+1,1]),n(i),[]);
        % dga (and possibly epsilon) are col vectors the same length as the
        % current dimension; copy this vector so it is the same size as
        % d2qdx2
        d2coef = repmat((epsilon.*dga{i}/2),[1,size(d2qdx2,2)]);
        % Apply ODE and permute in reverse order:
        dq = permute(reshape(dt*(-dfdx + d2coef.*d2qdx2),s),[dim+1, unrot]);
        
        % Apply BC's
        % This needs to be vectorized, it's really slow
        names = fields(bc);
        for bc_num=1:length(names)
            array = bc.(names{bc_num});
            nodes = array{1};
            stencil = array{2};
            weights = array{3};
            
            for node_index=1:size(nodes,2)
                total = zeros(n_eq,1);
                for stencil_index=1:size(stencil,2)
                    % find node by adding node coords to stencil coords
                    coords = nodes(:,node_index)+stencil(:,stencil_index);
                    temp(1) = {':'};
                    if dim == 1
                        temp(2) = {coords(1)};
                    elseif dim == 2
                        temp(2:3) = {coords(1);coords(2)};
                    elseif dim == 3
                        temp(2:4) = {coords(1);coords(2);coords(3)};
                    else
                        temp(2:dim+1) = num2cell(coords);
                    end
                    total = total + weights(:,:,stencil_index)*qc(temp{:});
                end
                coords = nodes(:,node_index);
                temp(1) = {':'};
                if dim == 1
                    temp(2) = {coords(1)};
                elseif dim == 2
                    temp(2:3) = {coords(1);coords(2)};
                elseif dim == 3
                    temp(2:4) = {coords(1);coords(2);coords(3)};
                else
                    temp(2:dim+1) = num2cell(coords);
                end
                dq(temp{:}) = dt*total;
            end
        end
        
        % Update the current state
        qc = qc + dq;
    end
    tc = tc + dt;
    iter = iter + 1;
end

end

function qnew = characteristic(amat, qc, d, dt)
eigvect = zeros(3,3,n);
eigval = zeros(3,3,n);
for i=1:n
    [eigvect(:,:,i),eigval(:,:,i)] = eig(amat(i));
end

% Characteristic variables
wc = zeros(3,n);
for i=1:n
    wc(:,i) = eigvect(:,:,i)\qc(:,i);
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
qnew = zeros(3,n);
for i=1:n
    qnew(:,i) = eigvect(:,:,i)*wn(:,i);
end
end


function d = central_fd(x,bc)
x = x(:);
dxp = circshift(x,-1)-x;
dxp(end) = nan; % dx to the right of points, last value is nan
dxm = circshift(dxp,1); % dx to the left of points, first value is nan

denom = dxp + dxm;
left = -(dxp./dxm)./denom;
center = ((dxp.^2 - dxm.^2)./(dxp.*dxm))./denom;
right = (dxm./dxp)./denom;

if strcmp(bc,'wall')
    lgc = logical(isnan(center));
    
    left(lgc) = 0; % Do not compute outer "ghost" nodes, will be set manually
    center(lgc) = 0;
    right(lgc) = 0;
elseif strcmp(bc,'open')
    error('Open not implemented');
elseif strcmp(bc,'periodic')
    dximg = (dxp(end-1)+dxp(1))/2; % imaginary dx between last point and first point of next period
    dxp(end) = dximg;
    dxm(1) = dximg;
    lgc = logical(isnan(center));
    
    denom(lgc) = dxp(lgc) + dxm(lgc);
    left(lgc) = -(dxp(lgc)./dxm(lgc))./denom(lgc);
    center(lgc) = ((dxp(lgc).^2 - dxm(lgc).^2)./(dxp(lgc).*dxm(lgc)))./denom(lgc);
    right(lgc) = (dxm(lgc)./dxp(lgc))./denom(lgc);

else
    error('Incorrect Boundary Condition');
end

d = (circshift(diag(right),1) + diag(center) + circshift(diag(left),-1))';

end

function d = central_fd_second(x,bc)
x = x(:);
dxp = circshift(x,-1)-x;
dxp(end) = nan; % dx to the right of points, last value is nan
dxm = circshift(dxp,1); % dx to the left of points, first value is nan

denom = dxp + dxm;
left = (2./dxm)./denom;
center = -2./(dxp.*dxm);
right = (2./dxp)./denom;

if strcmp(bc,'wall')
    lgc = logical(isnan(center));
    
    left(lgc) = 0; % Do not compute outer "ghost" nodes, will be set manually
    center(lgc) = 0;
    right(lgc) = 0;
elseif strcmp(bc,'open')
    error('Open not implemented');
elseif strcmp(bc,'periodic')
    dximg = (dxp(end-1)+dxp(1))/2; % imaginary dx between last point and first point of next period
    dxp(end) = dximg;
    dxm(1) = dximg;
    lgc = logical(isnan(center));
    
    denom(lgc) = dxp(lgc) + dxm(lgc);
    left(lgc) = (2./dxm(lgc))./denom(lgc);
    center(lgc) = -2./(dxp(lgc).*dxm(lgc));
    right(lgc) = (2./dxp(lgc))./denom(lgc);
else
    d = nan(length(dx));
    error('Incorrect Boundary Condition');
end
d = (circshift(diag(right),1) + diag(center) + circshift(diag(left),-1))';

end