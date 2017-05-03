function g = EulerFunY2D(q,~,gamma)
g = zeros(size(q));
q1 = q(1,:,:);
q2 = q(2,:,:);
q3 = q(3,:,:);
q4 = q(4,:,:);

v = sqrt(q2.^2 + q3.^2)./q1;
p = (gamma - 1)*(q4 - (1/2)*q1.*v.^2);

g(1,:,:) = q3;
g(2,:,:) = q3.*q2./q1;
g(3,:,:) = (q3.^2)./q1 + p;
g(4,:,:) = q3.*(q4 + p)./q1;
end