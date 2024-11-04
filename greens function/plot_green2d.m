% %% Plotting free space green's function in 2d
% %
% % 
% clf 
% clear
% close all

h = 0.2;
xs = -15:h:15;
[X,Y] = meshgrid(xs);
targ = [X(:).'; Y(:).'];
src = [0; 0];
gamma = -1; % fails if zero
beta = 3;


[rts,ejs] = find_roots(beta,gamma);

gf = green(src,targ,rts,ejs);
val = gf{1};
hess = gf{2};
gradlap = gf{3};
phi = gf{4};

val = reshape(val,size(X));
phi = reshape(phi,size(X));
hessxx = reshape(hess(:,:,1),size(X)); 
hessxy = reshape(hess(:,:,2),size(X)); 
hessyy = reshape(hess(:,:,3),size(X)); 
lapx = reshape(gradlap(:,:,1),size(X)); 
lapy = reshape(gradlap(:,:,2),size(X)); 

% figure(1);
% t = tiledlayout(1,3);
% title(t, ['Greens function, \beta = ',num2str(beta), ', \gamma = ',num2str(gamma)]);
% 
% nexttile
% surf = pcolor(X,Y,real(val));
% surf.EdgeColor = 'none';
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(val));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(val));
% surf.EdgeColor = 'none';
% colorbar
% title('abs')
% 
% figure(2);
% t = tiledlayout(1,3);
% title(t, ['G_\phi, \beta = ',num2str(beta), ', \gamma = ',num2str(gamma)]);
% 
% nexttile
% surf = pcolor(X,Y,real(phi));
% surf.EdgeColor = 'none';
% colorbar
% title('real')
% 
% nexttile
% surf = pcolor(X,Y,imag(phi));
% surf.EdgeColor = 'none';
% colorbar
% title('imag')
% 
% nexttile
% surf = pcolor(X,Y,abs(phi));
% surf.EdgeColor = 'none';
% colorbar
% title('abs') 


%% Second derivatives


figure(5);
t = tiledlayout(1,3);
title(t, 'G_{xx}');

nexttile
surf = pcolor(X,Y,real(hessxx)); % ;
surf.EdgeColor = 'none';
colorbar
title('real')

nexttile
surf = pcolor(X,Y,imag(hessxx));
surf.EdgeColor = 'none';
colorbar
title('abs')

nexttile
surf = pcolor(X,Y,abs(hessxx));
surf.EdgeColor = 'none';
colorbar
title('abs')



figure(6);
t = tiledlayout(1,3);
title(t, 'G_{xy}');

nexttile
surf = pcolor(X,Y,real(hessxy));
surf.EdgeColor = 'none';
colorbar
title('real')

nexttile
surf = pcolor(X,Y,imag(hessxy));
surf.EdgeColor = 'none';
colorbar
title('imag')

nexttile
surf = pcolor(X,Y,abs(hessxy));
surf.EdgeColor = 'none';
colorbar
title('abs')

figure(7)
t = tiledlayout(1,3);
title(t, 'G_{yy}');

nexttile
surf = pcolor(X,Y,real(hessyy));
surf.EdgeColor = 'none';
colorbar
title('real')

nexttile
surf = pcolor(X,Y,imag(hessyy));
surf.EdgeColor = 'none';
colorbar
title('imag')

nexttile
surf = pcolor(X,Y,abs(hessyy));
surf.EdgeColor = 'none';
colorbar
title('abs')


%% Third derivatives 

lap = hessxx + hessyy;

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