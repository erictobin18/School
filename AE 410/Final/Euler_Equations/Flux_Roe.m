function flux = Flux_Roe(q_l,q_r, ~, g)
if length(q_l)~=length(q_r)
    error('Tried to find fluxes between unequal length vectors')
end
n = size(q_l,1);
% find averaged values
uhat = (q_l(:,1).^(-1/2).*q_l(:,2) + q_r(:,1).^(-1/2).*q_r(:,2))./(sqrt(q_l(:,1))+sqrt(q_r(:,1)));
hhat = (sqrt(q_l(:,1)).*enthalpy(q_l,g) + sqrt(q_r(:,1)).*enthalpy(q_r,g))./(sqrt(q_l(:,1))+sqrt(q_r(:,1)));
chat = (g-1)*(hhat - (uhat.^2)/2);
rhohat = (q_l(:,1).*q_r(:,1)).^(1/2);
eigvals = [ uhat, uhat - chat, uhat + chat ];
epsilon = ones(n,1)*.001;
L2 = logical(abs(eigvals(:,2))<epsilon); % Guarantee abs(lambda) > 0
L3 = logical(abs(eigvals(:,3))<epsilon);
eigvals(L2,2) = (1/2)*((eigvals(L2,2).^2)./epsilon(L2,1)+epsilon(L2,1));
eigvals(L3,3) = (1/2)*((eigvals(L3,3).^2)./epsilon(L3,1)+epsilon(L3,1));

one = ones(n,1);
eigvect1 = [ one, uhat, (uhat.^2)/2];
eigvect2 = [ one, uhat-chat, hhat-uhat.*chat];
eigvect3 = [ one, uhat+chat, hhat+uhat.*chat];
drho = q_r(:,1) - q_l(:,1);
dp = ideal_p(q_r,g) - ideal_p(q_l,g);
du = q_r(:,2)./q_r(:,1) - q_l(:,2)./q_l(:,1);
w = [ drho - dp./chat.^2, (dp - chat.*rhohat.*du)./(2*chat.^2),(dp + chat.*rhohat.*du)./(2*chat.^2) ];

avgterm = (1/2)*[ q_l(:,2) + q_r(:,2), (q_l(:,2).^2)./q_l(:,1) + ideal_p(q_l,g) + (q_r(:,2).^2)./q_r(:,1) + ideal_p(q_r,g), q_l(:,2).*(q_l(:,3)+ideal_p(q_l,g))./q_l(:,1) + q_r(:,2).*(q_r(:,3)+ideal_p(q_r,g))./q_r(:,1) ];
term1 = (1/2)*(repmat(abs(eigvals(:,1)),1,3).*eigvect1.*repmat(w(:,1),1,3));
term2 = (1/2)*(repmat(abs(eigvals(:,2)),1,3).*eigvect2.*repmat(w(:,2),1,3));
term3 = (1/2)*(repmat(abs(eigvals(:,3)),1,3).*eigvect3.*repmat(w(:,3),1,3));
flux = avgterm - term1 - term2 - term3;
end
function out = enthalpy(y,g)
out = y(:,3)./y(:,1) - (1/2)*(y(:,2)./y(:,1)).^2 + ideal_p(y,g)./y(:,1);
end
function out = ideal_p(y,g) % solve equation of state for p
out = (g - 1)*(y(:,3) - (1/2)*(y(:,2).^2)./y(:,1));
end