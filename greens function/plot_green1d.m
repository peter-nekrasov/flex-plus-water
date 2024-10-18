%% Plotting green's function in 1d
%

beta = 2;
gamma = -0.5;
h = 0.01;
xs = -20:h:20;
[val, grad, hess] = green(xs,0*xs,beta,gamma,false);

grad = grad(:,:,1);
hess = hess(:,:,1);

figure(8)
plot(xs,real(hess),xs, imag(hess))
legend('real','imaginary')
title(['greens function, \beta = ', num2str(beta),', \gamma = ', num2str(gamma)])

val = gradient(grad, h);

figure(9)
plot(xs,real(val),xs, imag(val))

