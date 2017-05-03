clear all; clc;

% run lecture7_1.m
lecture7_1

% compute continuous eigenvalues (positive imaginary components only)
for m = 0:(N-1)/2
  lambda_c(m+1,1) = 1i*(2*pi/L)*m;
end

% discrete eigenvalues (positive imaginary components only)
lambda = lambda(1:(N-1)/2+1);


% plot eigenvalues
figure
subplot(1,2,1);

plot(real(lambda_c),imag(lambda_c),'o'); hold on;
plot(real(lambda),  imag(lambda)  ,'x'); 
xlim([-max(max(abs(real(lambda))),1) max(max(abs(real(lambda))),1)]);

xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
legend('Continuous eigenvalues','Discrete eigenvalues',...
       'location','southoutside');

subplot(1,2,2);

plot(imag(lambda_c),'o'); hold on;
plot(imag(lambda),  'x')

xlabel('m'); ylabel('Im(\lambda_m)');
legend('Continuous eigenvalues','Discrete eigenvalues',...
       'location','southoutside');
