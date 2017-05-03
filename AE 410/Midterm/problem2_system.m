function dqdt = problem2_system(t,q,D,D2);
N = size(D,1);

u_1 = q(1:N);
u_2 = q(N+1:2*N);
u_3 = q(2*N+1:3*N);

du_1dt =  0*D*u_1 - 1*D*u_2 -  0*D*u_1 + D2*u_1; 
du_2dt = .5*D*u_1 - 1*D*u_2 - .5*D*u_3 + D2*u_2;
du_3dt =  3*D*u_1 - 3*D*u_2 -  1*D*u_3 + D2*u_3;

dqdt = [du_1dt; du_2dt; du_3dt];
