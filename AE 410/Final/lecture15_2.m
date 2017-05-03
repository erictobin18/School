clear all;
clc

% check if flux is homogeneous of degree 1 -----------------------------------
gamma = 1.4;

F = @(q) [q(2);...
          q(2)^2/q(1) + (gamma-1)*(q(3)-0.5*q(2)^2/q(1));...
          q(3) + (gamma-1)*(q(3) - 0.5*q(2)^2/q(1))*(q(2)/q(1))];

% if true, the F(\lambda*q) = lambda*F(q) for arbitrary lamda
q = rand(3,1);
lambda = rand;
F(lambda*q) - lambda*F(q)

% check if F=A(q)*q ----------------------------------------------------------
syms gamma;
q = sym('q_%d', [3 1]);

F = [q(2);...
     q(2)^2/q(1) + (gamma-1)*(q(3)-0.5*q(2)^2/q(1));...
     q(3) + (gamma-1)*(q(3) - 0.5*q(2)^2/q(1))*(q(2)/q(1))];

A = simplify(jacobian(F,q),'Steps',200);

simplify(F - A*q,'Steps',200)
