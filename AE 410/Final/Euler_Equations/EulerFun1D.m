function f = EulerFun1D(q, ~,gamma)
f = zeros(size(q));
q1 = q(1,:);
q2 = q(2,:);
q3 = q(3,:);

u = q2./q1;
p = (gamma - 1)*(q3 - (1/2)*q2.*u);

f(1,:) = q2;
f(2,:) = q2.*u + p;
f(3,:) = u.*(q3 + p);
end