function dqdt = problem2_system_upwind(t,q,Ap,Am,Dp,Dm);
N = size(Dp,1);

u_1 = q(1:N);
u_2 = q(N+1:2*N);
u_3 = q(2*N+1:3*N);



du_1dt = - Ap(1,1)*Dp*u_1 - Ap(1,2)*Dp*u_2 - Ap(1,3)*Dp*u_3 - Am(1,1)*Dm*u_1 - Am(1,2)*Dm*u_2 - Am(1,3)*Dm*u_3; 
du_2dt = - Ap(2,1)*Dp*u_1 - Ap(2,2)*Dp*u_2 - Ap(2,3)*Dp*u_3 - Am(2,1)*Dm*u_1 - Am(2,2)*Dm*u_2 - Am(2,3)*Dm*u_3;
du_3dt = - Ap(3,1)*Dp*u_1 - Ap(3,2)*Dp*u_2 - Ap(3,3)*Dp*u_3 - Am(3,1)*Dm*u_1 - Am(3,2)*Dm*u_2 - Am(3,3)*Dm*u_3;
dqdt = [du_1dt; du_2dt; du_3dt];
