% Plotting green's function
%

beta = 0.5;
gamma = -1;
xs = 0:0.01:50*beta^(-1/4);

plot(xs,real(green(xs,beta,gamma)),xs, imag(green(xs,beta,gamma)))
legend('real','imaginary')
title(['greens function, \beta = ', num2str(beta),', \gamma = ', num2str(gamma)])