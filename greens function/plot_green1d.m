%% Plotting green's function in 1d
%

beta = 2;
gamma = -0.5;
h = 0.01;
xs = -20:h:20;

[rts,ejs] = find_roots(beta,gamma);

gf = green([0; 0],[xs; 0*xs],rts,ejs);

val = gf{1};
hess = gf{2};
third = gf{3};
phi = gf{4};

tiledlayout(1,4)
nexttile
plot(xs,real(val),xs, imag(val))
legend('real','imaginary')
title('G')

nexttile
plot(xs,real(hess(:,:,1)),xs, imag(hess(:,:,1)))
legend('real','imaginary')
title('\partial_{xx} G')

nexttile
plot(xs,real(third(:,:,1)),xs, imag(third(:,:,1)))
legend('real','imaginary')
title('\partial_{xxx} G')

nexttile
plot(xs,real(phi),xs, imag(phi))
legend('real','imaginary')
title('G_\phi')
