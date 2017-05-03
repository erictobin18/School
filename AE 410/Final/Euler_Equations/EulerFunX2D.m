function f = EulerFunX2D(q,~,gamma)
f = zeros(size(q));
q1 = q(1,:,:);
q2 = q(2,:,:);
q3 = q(3,:,:);
q4 = q(4,:,:);

v = sqrt(q2.^2 + q3.^2)./q1;
p = (gamma - 1)*(q4 - (1/2)*q1.*v.^2);

f(1,:,:) = q2;
f(2,:,:) = (q2.^2)./q1 + p;
f(3,:,:) = q2.*q3./q1;
f(4,:,:) = q2.*(q4 + p)./q1;
end