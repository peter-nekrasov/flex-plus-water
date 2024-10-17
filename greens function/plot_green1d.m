%% Plotting green's function in 1d
%

beta = 2;
gamma = -0.5;
xs = -20:0.01:20;
[val, grad] = green(xs,0*xs,beta,gamma,false);

figure(8)
plot(xs,real(val),xs, imag(val))
legend('real','imaginary')
title(['greens function, \beta = ', num2str(beta),', \gamma = ', num2str(gamma)])

disp(val(10))

figure(9)
plot(xs,real(grad(:,:,1)),xs, imag(grad(:,:,1)))

