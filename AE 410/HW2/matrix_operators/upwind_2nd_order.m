function D = upwind_2nd_order(x,PeriodicFlag)

  % Returns the second-order, upwind finite difference 
  % approximation of the first derivative

Nx = max(size(x));
dx = x(2) - x(1);

V0 = ones([1,Nx]);
V1 = ones([1,Nx-1]);
V2 = ones([1,Nx-2]);
V3 = ones([1,Nx-3]);

D = (diag(V2,-2) - 4*diag(V1,-1) + 3*diag(V0,0))/(2*dx);

if(PeriodicFlag == 0)
  D(1,:)  = [-1 1 0*V2]/(dx);
  D(2,:)  = [-1 1 0*V2]/(dx);
elseif(PeriodicFlag == 1) 
  D(1,:)    = [3 0*V3 1 -4]/(2*dx);
  D(2,:)    = [-4 3 0*V3 1]/(2*dx);
else
  error('Incorrect PeriodicFlag');
end
