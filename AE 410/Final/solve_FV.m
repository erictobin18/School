function [tsol, qsol] = solve_FV(grid, t, f, q0, bc) % implementation space

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
clear('temp');

qc = q0;
qsol = zeros([n_eq, n', nt]); % save each variable at every grid point at every timestep

iter = 1; % iteration number

% flux from previous timestep stored with same shape as q
while tc < tf + dt
    if mod(iter,10)==0
        disp(strcat('Iteration ',num2str(iter)));
    end
    qsol(:,g_shape{:},iter) = qc; % save qc
    tsol(iter) = tc; % save time of current iteration
    
    
    
    for i=1:dim % Lie Splitting: Solve each dimension independently first-order
        % transpose the current values so the dimension of interest is
        % along the cols:
        rot = circshift((1:dim)',1-i)';
        unrot = circshift((1:dim)',i-1)';
        qtemp = permute(qc,[rot + 1, 1]);
                          % col
        f_r = f{i}(qtemp,circshift(qtemp,-1),t); % array of fluxes to the "right" of each point
%         f_r = zeros(size(qtemp,1),3);
%         cs = circshift(qtemp,-1);
% %         tmp = zeros(size(qtemp,1));
%         for j=1:size(qtemp,1)
%             f_r(j,:) = AppRiemannSol(qtemp(j,:)',cs(j,:)')';
%         end
%         f_r = f_ri;
        f_l = circshift(f_r,1); % array of fluxes to the left of each point
        
        s = size(f_r);
        dg = repmat(dga{i},[1,s(2:end)]);
        dq = permute((f_l - f_r)*dt./dg,[dim+1,unrot]);
        
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
                    total = total + weights(:,:,stencil_index)*dq(temp{:});
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
        % store the fluxes for the next timestep
        
    end
    tc = tc + dt;
    iter = iter + 1;
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
