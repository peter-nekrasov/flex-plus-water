% Plotting free space green's function in 2d
%

h = 0.1;
xs = -10:h:10;
[X,Y] = meshgrid(xs);
gamma = -1; % fails if zero
beta = 10;
R = sqrt(X.^2 + Y.^2);
gval = green(R,beta,gamma,false);

t = tiledlayout(1,3);
title(t, ['\beta = ',num2str(beta), ', \gamma = ',num2str(gamma)]);

nexttile
h = pcolor(X,Y,real(gval));
h.EdgeColor = 'none';
colorbar
title('real')

nexttile
h = pcolor(X,Y,imag(gval));
h.EdgeColor = 'none';
colorbar
title('imag')

nexttile
h = pcolor(X,Y,abs(gval));
h.EdgeColor = 'none';
colorbar
title('abs')