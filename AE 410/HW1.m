function HW1

dx = 10.^(-10:.1:2);
xi = .1*ones(1,length(dx));
xim1 = xi - dx;
xim2 = xi - 2*dx;

ui = sin(xi);
uim1 = sin(xim1);
uim2 = sin(xim2);

a1 = 1./dx; b1 = -1./dx;
err1 = abs(cos(xi) - a1.*ui - b1.*uim1);

a2 = 3./(2*dx); b2 = -2./(dx); c2 = 1./(2*dx);
err2 = abs(cos(xi) - a2.*ui - b2.*uim1 - c2.*uim2);

loglog(dx,err1,'-k',dx,err2,'-.k');
grid on;
axis on;
xlabel('Delta X');
ylabel('Error');
title('Scheme Error vs. Delta X');
legend('show');
end

