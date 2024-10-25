% %% Plotting free space green's function in 2d
% %
% % 
% clf 
% clear
% close all

h = 0.05;
xs = -10:h:10;
[X,Y] = meshgrid(xs);
gamma = -1; % fails if zero
beta = 10;

gf = green(X,Y,beta,gamma,false);

figure(1);
t = tiledlayout(1,3);
title(t, ['\beta = ',num2str(beta), ', \gamma = ',num2str(gamma)]);

nexttile
surf = pcolor(X,Y,real(gf{1}));
surf.EdgeColor = 'none';
colorbar
title('real')

nexttile
surf = pcolor(X,Y,imag(gf{1}));
surf.EdgeColor = 'none';
colorbar
title('imag')

nexttile
surf = pcolor(X,Y,abs(gf{1}));
surf.EdgeColor = 'none';
colorbar
title('abs')


%% Second derivatives


figure(5);
t = tiledlayout(1,3);
title(t, 'G_{xx}');

nexttile
surf = pcolor(X,Y,real(gf{2}));
surf.EdgeColor = 'none';
clim([-1.2 0.3])
colorbar
title('real')

nexttile
surf = pcolor(X,Y,imag(gf{2}));
surf.EdgeColor = 'none';
colorbar
title('imag')

nexttile
surf = pcolor(X,Y,abs(gf{2}));
surf.EdgeColor = 'none';
colorbar
title('abs')


figure(6);
t = tiledlayout(1,3);
title(t, 'G_{xy}');

nexttile
surf = pcolor(X,Y,real(gf{3}));
surf.EdgeColor = 'none';
colorbar
title('real')

nexttile
surf = pcolor(X,Y,imag(gf{3}));
surf.EdgeColor = 'none';
colorbar
title('imag')

nexttile
surf = pcolor(X,Y,abs(gf{3}));
surf.EdgeColor = 'none';
colorbar
title('abs')

figure(7)
t = tiledlayout(1,3);
title(t, 'G_{yy}');

nexttile
surf = pcolor(X,Y,real(gf{4}));
surf.EdgeColor = 'none';
colorbar
title('real')

nexttile
surf = pcolor(X,Y,imag(gf{4}));
surf.EdgeColor = 'none';
colorbar
title('imag')

nexttile
surf = pcolor(X,Y,abs(gf{4}));
surf.EdgeColor = 'none';
colorbar
title('abs')


%% Third derivatives 

lap = gf{2} + gf{4};
lapx = gf{5};
lapy = gf{6};

figure(11);
t = tiledlayout(1,3);
title(t, '\Delta G');

nexttile
surf = pcolor(X,Y,real(lap));
surf.EdgeColor = 'none';
colorbar
title('real')

nexttile
surf = pcolor(X,Y,imag(lap));
surf.EdgeColor = 'none';
colorbar
title('imag')

nexttile
surf = pcolor(X,Y,abs(lap));
surf.EdgeColor = 'none';
colorbar
title('abs')


figure(12);
t = tiledlayout(1,3);
title(t, '\partial_x \Delta G');

nexttile
surf = pcolor(X,Y,real(lapx));
surf.EdgeColor = 'none';
colorbar
title('real')

nexttile
surf = pcolor(X,Y,imag(lapx));
surf.EdgeColor = 'none';
colorbar
title('imag')

nexttile
surf = pcolor(X,Y,abs(lapx));
surf.EdgeColor = 'none';
colorbar
title('abs')


figure(13);
t = tiledlayout(1,3);
title(t, '\partial_y \Delta G');

nexttile
surf = pcolor(X,Y,real(lapy));
surf.EdgeColor = 'none';
colorbar
title('real')

nexttile
surf = pcolor(X,Y,imag(lapy));
surf.EdgeColor = 'none';
colorbar
title('imag')

nexttile
surf = pcolor(X,Y,abs(lapy));
surf.EdgeColor = 'none';
colorbar
title('abs')

[lapx,lapy] = gradient(lap,h);

figure(14);
t = tiledlayout(1,3);
title(t, 'numerical \partial_x \Delta G');

nexttile
surf = pcolor(X,Y,real(lapx));
surf.EdgeColor = 'none';
colorbar
title('real')

nexttile
surf = pcolor(X,Y,imag(lapx));
surf.EdgeColor = 'none';
colorbar
title('imag')

nexttile
surf = pcolor(X,Y,abs(lapx));
surf.EdgeColor = 'none';
colorbar
title('abs')


figure(15);
t = tiledlayout(1,3);
title(t, 'numerical \partial_y \Delta G');

nexttile
surf = pcolor(X,Y,real(lapy));
surf.EdgeColor = 'none';
colorbar
title('real')

nexttile
surf = pcolor(X,Y,imag(lapy));
surf.EdgeColor = 'none';
colorbar
title('imag')

nexttile
surf = pcolor(X,Y,abs(lapy));
surf.EdgeColor = 'none';
colorbar
title('abs')