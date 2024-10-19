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

% [gval,grad] = green(X,Y,beta,gamma,false);
% 
% figure(1);
% t = tiledlayout(1,3);
% title(t, ['\beta = ',num2str(beta), ', \gamma = ',num2str(gamma)]);
% 
% nexttile
% surf = pcolor(X,Y,real(gval));
% surf.EdgeColor = 'none';
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(gval));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(gval));
% surf.EdgeColor = 'none';
% colorbar
% title('abs')
% 
% %% First derivatives
% 
% figure(2);
% t = tiledlayout(1,3);
% title(t, ['\beta = ',num2str(beta), ', \gamma = ',num2str(gamma)]);
% 
% nexttile
% surf = pcolor(X,Y,real(grad(:,:,1)));
% surf.EdgeColor = 'none';
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(grad(:,:,1)));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(grad(:,:,1)));
% surf.EdgeColor = 'none';
% colorbar
% title('abs')
% 
% figure(3);
% t = tiledlayout(1,3);
% title(t, ['\beta = ',num2str(beta), ', \gamma = ',num2str(gamma)]);
% 
% nexttile
% surf = pcolor(X,Y,real(grad(:,:,2)));
% surf.EdgeColor = 'none';
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(grad(:,:,2)));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(grad(:,:,2)));
% surf.EdgeColor = 'none';
% colorbar
% title('abs')
% 
% %% Checking against numerical gradient
% 
% [FX,FY] = gradient(gval,h);
% 
% figure(4);
% t = tiledlayout(1,3);
% title(t, ['\beta = ',num2str(beta), ', \gamma = ',num2str(gamma)]);
% 
% nexttile
% surf = pcolor(X,Y,real(FX));
% surf.EdgeColor = 'none';
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(FX));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(FX));
% surf.EdgeColor = 'none';
% colorbar
% title('abs')


% Second derivatives


% [gval,grad,hess] = green(X,Y,beta,gamma,false);
% 
% 
% figure(5);
% t = tiledlayout(1,3);
% title(t, 'G_{xx}');
% 
% nexttile
% surf = pcolor(X,Y,real(hess(:,:,1)));
% surf.EdgeColor = 'none';
% clim([-1.2 0.3])
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(hess(:,:,1)));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(hess(:,:,1)));
% surf.EdgeColor = 'none';
% colorbar
% title('abs')
% 
% 
% figure(6);
% t = tiledlayout(1,3);
% title(t, 'G_{xy}');
% 
% nexttile
% surf = pcolor(X,Y,real(hess(:,:,2)));
% surf.EdgeColor = 'none';
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(hess(:,:,2)));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(hess(:,:,2)));
% surf.EdgeColor = 'none';
% colorbar
% title('abs')
% 
% figure(7)
% t = tiledlayout(1,3);
% title(t, 'G_{yy}');
% 
% nexttile
% surf = pcolor(X,Y,real(hess(:,:,3)));
% surf.EdgeColor = 'none';
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(hess(:,:,3)));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(hess(:,:,3)));
% surf.EdgeColor = 'none';
% colorbar
% title('abs')
% 
% % Calculating second derivatives numerically 
% 
% [FXX,FXY] = gradient(grad(:,:,1),h);
% [~,FYY] = gradient(grad(:,:,2),h);
% 
% 
% figure(8);
% t = tiledlayout(1,3);
% title(t, 'Numerical G_{xx}');
% 
% nexttile
% surf = pcolor(X,Y,real(FXX));
% surf.EdgeColor = 'none';
% clim([-1.2 0.3])
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(FXX));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(FXX));
% surf.EdgeColor = 'none';
% colorbar
% title('abs')
% 
% figure(9);
% t = tiledlayout(1,3);
% title(t, 'Numerical G_{xy}');
% 
% nexttile
% surf = pcolor(X,Y,real(FXY));
% surf.EdgeColor = 'none';
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(FXY));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(FXY));
% surf.EdgeColor = 'none';
% colorbar
% title('abs')
% 
% %% Checking Laplacian
% 
% lap = FXX + FYY;
% 
% figure(10);
% t = tiledlayout(1,3);
% title(t, 'numerical \Delta G');
% 
% nexttile
% surf = pcolor(X,Y,real(lap));
% surf.EdgeColor = 'none';
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(lap));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(lap));
% surf.EdgeColor = 'none';
% colorbar
% title('abs')
% 
% lap = hess(:,:,1) + hess(:,:,3);
% 
% figure(11);
% t = tiledlayout(1,3);
% title(t, '\Delta G');
% 
% nexttile
% surf = pcolor(X,Y,real(lap));
% surf.EdgeColor = 'none';
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(lap));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(lap));
% surf.EdgeColor = 'none';
% colorbar
% title('abs')

%% Third derivatives 

[~,~,hess,gradlap] = green(X,Y,beta,gamma,false);

lapx = gradlap(:,:,1);
lapy = gradlap(:,:,2);

lap = hess(:,:,1) + hess(:,:,3);

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