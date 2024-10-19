%% Plotting green's function in 1d
%

beta = 2;
gamma = -0.5;
h = 0.01;
xs = -20:h:20;
[val, grad, hess,third] = green(xs,0*xs,beta,gamma);

grad = grad(:,:,1);
hess = hess(:,:,1);
third = third(:,:,1);

figure(1)
plot(xs,real(val),xs, imag(val))
legend('real','imaginary')
title(['greens function, \beta = ', num2str(beta),', \gamma = ', num2str(gamma)])



figure(2)
plot(xs,real(grad),xs, imag(grad))
legend('real','imaginary')
title('\partial_x G')


figure(3)
plot(xs,real(hess),xs, imag(hess))
legend('real','imaginary')
title('\partial_{xx} G')

figure(4)
plot(xs,real(third),xs, imag(third))
legend('real','imaginary')
title('\partial_{xxx} G')

