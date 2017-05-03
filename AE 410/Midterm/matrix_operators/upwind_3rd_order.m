function D = upwind_3rd_order(x,PeriodicFlag)

  % Returns the third-order, upwind finite difference 
  % approximation of the first derivative

Nx = max(size(x));
dx = x(2) - x(1);

V0 = ones([1,Nx]);
V1 = ones([1,Nx-1]);
V2 = ones([1,Nx-2]);
V3 = ones([1,Nx-3]);
V4 = ones([1,Nx-4]);

D = (diag(V2,-2) - 6*diag(V1,-1) + 3*diag(V0,0) + 2*diag(V1,+1))/(6*dx);

if(PeriodicFlag == 0)
  D(1,:)  = [-1 1 0*V2]/(dx);
  D(2,:)  = [-1 1 0*V2]/(dx);
  D(Nx,:) = [0*V2 -1 1]/(dx);
elseif(PeriodicFlag == 1) 
  D(1,:)  = [3 2 0*V4 1 -6]/(6*dx);
  D(2,:)  = [-6 3 2 0*V4 1]/(6*dx);
  D(Nx,:) = [2 0*V4 1 -6 3]/(6*dx);
else
  error('Incorrect PeriodicFlag');
end
