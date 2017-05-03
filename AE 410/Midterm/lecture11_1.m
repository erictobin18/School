clear all;
clc

% declare symbolic variables (g = gamma)
syms rho u p g c;

% assume gamma, pressure and density are positive
assume([g p rho], 'positive');

% speed of sound
c = sqrt((g*p)/rho);

A = [u    rho   0       ; 
     0    u     rho^-1  ;
     0    g*p   u       ];

X = [1  1/c^2           1/c^2             ;
     0  1/sqrt(g*p*rho) -1/sqrt(g*p*rho)  ;
     0  1               1                 ];

% Compare with Eq. 10.10a -- 10.10c
Lambda = inv(X)*A*X;
Lambda = simplify(Lambda,'steps',200);
fprintf('Lambda = \n')
pretty(Lambda)

% Compare with Eqs. 10.11
w = inv(X)*[rho; u; p];
w = simplify(w,'steps',200);
fprintf('w = \n')
pretty(w)
