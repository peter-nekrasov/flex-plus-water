% Plotting green's function in 1d
%

beta = 5;
gamma = -1;
xs = 0:0.01:50;
val = green(xs,beta,gamma);

figure(8)
plot(xs,real(val),xs, imag(val))
legend('real','imaginary')
title(['greens function, \beta = ', num2str(beta),', \gamma = ', num2str(gamma)])

disp(val(10))