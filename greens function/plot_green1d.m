% Plotting green's function in 1d
%

beta = 2;
gamma = -0.5;
xs = 0:0.1:20;
val = green(xs,beta,gamma,false);

figure(8)
plot(xs,real(val),xs, imag(val))
legend('real','imaginary')
title(['greens function, \beta = ', num2str(beta),', \gamma = ', num2str(gamma)])

disp(val(10))