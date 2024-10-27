%% Plotting green's function and comparing to flexural wave greens function
%

k = 2;
beta = k^4;
gamma = 0;
h = 0.01;
xs = -20:h:20;
gf = green(xs,0*xs,beta,gamma);
val = gf{1} / (2*pi);

figure(1)
tiledlayout(1,2)
nexttile
plot(xs,real(val),xs, imag(val))
legend('real','imaginary')
title(['greens function, \beta = ', num2str(beta),', \gamma = ', num2str(gamma)])

gf2 = 1/(2*k^2)*(1i/4*besselh(0,k*abs(xs)) - 1/(2*pi)*besselk(0,k*abs(xs)));

nexttile
plot(xs,real(val),xs, imag(val))
legend('real','imaginary')
title(['greens function, \beta = ', num2str(beta),', \gamma = ', num2str(gamma)])
