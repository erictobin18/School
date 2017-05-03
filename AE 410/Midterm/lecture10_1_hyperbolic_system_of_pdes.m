function dqdt = lecture10_1_hyperbolic_system_of_pdes(t,q,D);
N = size(D,1);

u_1 = q(1:N);
u_2 = q(N+1:2*N);
u_3 = q(2*N+1:3*N);

du_1dt = -2*D*u_1 - 4*D*u_2          ; 
du_2dt =            2*D*u_2 - 1*D*u_3;
du_3dt =          - 1*D*u_2 - 2*D*u_3;

dqdt = [du_1dt; du_2dt; du_3dt];
