function D = downwind_1st_order(x,PeriodicFlag)

  % Returns the first-order, downwind finite difference 
  % approximation of the first derivative 

Nx = max(size(x));
dx = x(2) - x(1);

V0 = ones([1,Nx]);
V1 = ones([1,Nx-1]);

D = (1/dx)*(-diag(V0,0) + diag(V1,+1));

if(PeriodicFlag == 0)
  D(end,end) = 1/dx;
  D(end,end-1) = -1/dx;
elseif(PeriodicFlag == 1) 
  D(end,1) = 1/dx;
else
  error('Incorrect PeriodicFlag');
end
