%% Plotting green's function in 1d
%

beta = 2;
gamma = -0.5;
h = 0.01;
xs = -10:h:10;

[rts,ejs] = find_roots(beta,gamma);

gf = green([0; 0],[xs; 0*xs],rts,ejs);

val = gf{1};
hess = gf{2};
third = gf{3};
phi = gf{4};

hessxx = hess(:,:,1);
hessyy = hess(:,:,3);
thirdx = third(:,:,1);

tiledlayout(1,4)
nexttile
plot(xs,real(val),xs, imag(val))
legend('real','imaginary')
title('G')

nexttile
plot(xs,real(hess(:,:,1)),xs, imag(hess(:,:,1)))
legend('real','imaginary')
title('\partial_{xx} G')
hold on

nexttile
plot(xs,real(third(:,:,1)),xs, imag(third(:,:,1)))
legend('real','imaginary')
title('\partial_{xxx} G')

nexttile
plot(xs,real(phi),xs, imag(phi))
legend('real','imaginary')
title('G_\phi')


%% Checking derivatives using finite difference

v = [-1/560	8/315 -1/5 8/5 -205/72 8/5 -1/5	8/315 -1/560] / h^2;
subval = val(650-4:650+4);
err1 = sum(v.*subval.') - hessxx(650)

v = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280] / h;
subval = hessxx(650-4:650+4) + hessyy(650-4:650+4);
err2 = sum(v.*subval.') - thirdx(650)
