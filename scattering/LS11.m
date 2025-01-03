%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for plane wave 
% scattering of flexural-gravity waves with a system of channels
%
%
%%%%%

clear 
close all
addpath(genpath('..'))

h = 0.5;
L = 70;

xs = -L/2:h:L/2;
xl = -L:h:L;

ys = -L/2:h:L/2;
yl = -L:h:L;

[~,~,coefs,Hs] = get_Y(-35,35,-5,65,h,-0.1,200);

[X,Y] = meshgrid(xs,ys);
[XL,YL] = meshgrid(xl,yl);


H = Hs{1};
Hx = Hs{2};
Hy = Hs{3};
Hxx = Hs{4};
Hxy = Hs{5};
Hyy = Hs{6};

figure(1);
tiledlayout(2,3);

nexttile
s = pcolor(X,Y,H);
s.EdgeColor = 'None';
colorbar
title('H')
drawnow

nexttile
s = pcolor(X,Y,Hx);
s.EdgeColor = 'None';
colorbar
title('H_x')
drawnow

nexttile
s = pcolor(X,Y,Hy);
s.EdgeColor = 'None';
colorbar
title('H_y')
drawnow


nexttile
s = pcolor(X,Y,Hxx);
s.EdgeColor = 'None';
colorbar
title('H_{xx}')
drawnow

nexttile
s = pcolor(X,Y,Hxy);
s.EdgeColor = 'None';
colorbar
title('H_{xy}')
drawnow

nexttile
s = pcolor(X,Y,Hyy);
s.EdgeColor = 'None';
colorbar
title('H_{yy}')
drawnow


a0 = coefs{1}; 
b0 = coefs{3}; 
g0 = coefs{5}; 

% Finding positive real roots
[rts,ejs] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0));
ejs = ejs/a0;


% RHS (Incident field)
k1 = k*cos(pi/3);
k2 = k*sin(pi/3);
phiinc = exp(1i*k1*X+1i*k2*Y);
phininc = k*exp(1i*k1*X+1i*k2*Y);
[rhs_vec, rhsp] = get_rhs_vec(coefs,k1,k2,phiinc);

figure(1);
tiledlayout(1,4);

nexttile
s = pcolor(X,Y,H);
s.EdgeColor = 'None';
colorbar
title('H')
drawnow

nexttile
s = pcolor(X,Y,(coefs{1} + coefs{2}));
s.EdgeColor = 'None';
colorbar
title('\alpha')
drawnow

nexttile
s = pcolor(X,Y,(coefs{2} + coefs{3}));
s.EdgeColor = 'None';
colorbar
title('\beta')
drawnow


nexttile
s = pcolor(X,Y,real(rhsp));
s.EdgeColor = 'None';
colorbar
title('rhs')
drawnow

% Constructing integral operators
src = [xl(ceil(end/2)); yl(ceil(end/2))];
targ = [XL(:).'; YL(:).'];

[inds,corrs] = get_correct(h,a0);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs), h, inds,corrs);
ind = find((XL == src(1)) & (YL == src(2)));
sz = size(XL);
kerns = gen_fft_kerns(kerns,sz,ind);
ekerns = {kerns{1}, kerns{4}};

% Solve with GMRES
start = tic;
mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,100,1e-9,200);
mu = reshape(mu, size(X));
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

%inds = find(ismember(Xeval, X).*ismember(Yeval, Y));
%mu_aug = Xeval*0;
%mu_aug(inds) = mu;

[phi, phi_n] = sol_eval_fft(mu,ekerns);

phiinc = exp(1i*k1*X+1i*k2*Y);
phininc = k*exp(1i*k1*X+1i*k2*Y);

phi_tot = phi + phiinc;
phi_n_tot = phi_n + phininc;

load gong.mat
sound(y)

%%


figure(2);
pc = pcolor(X,Y,abs(mu));
pc.EdgeColor = 'none';
title('Abs(\mu)')
colorbar

figure(3);
tiledlayout(2,2)

nexttile
pc = pcolor(X,Y,real(phi_tot.'));
pc.EdgeColor = 'none';
title('Re(\phi)')
colorbar

nexttile
pc = pcolor(X,Y,abs(phi_tot.'));
pc.EdgeColor = 'none';
title('|\phi|')
colorbar


nexttile
pc = pcolor(X,Y,real(phi_n_tot.'));
pc.EdgeColor = 'none';
title('real(\phi_n)')
clim([-4 4])
colorbar

nexttile
pc = pcolor(X,Y,abs(phi_n_tot.'));
pc.EdgeColor = 'none';
title('|\phi_n|')
colorbar
       
% Calculate error with finite difference
err = get_fin_diff_err(X,Y,mu,phi_n_tot,phi_tot,h,coefs,0.5,1.5)


return

%% Figure generation for Jeremy

f = figure(4);

t = tiledlayout(1,3,'TileSpacing','tight'); 


% ax1 = nexttile;
% s = pcolor(X1,Y1,H);
% s.EdgeColor = 'None';
% colormap(ax1,c)
% clim([min(H(:)) max(H(:))])
% colorbar
% title('Thickness (m)','FontWeight','normal')


nexttile
pc = pcolor(X,Y,real(H.'));
pc.EdgeColor = 'none';
title('H')
colorbar
axis square

nexttile
pc = pcolor(X,Y,real(phi_n_tot.'));
pc.EdgeColor = 'none';
title('re(\phi_n)')
clim([-4 4])
colorbar
axis square


nexttile
pc = pcolor(X,Y,abs(phi_n_tot.'));
pc.EdgeColor = 'none';
title('|\phi_n|')
colorbar
axis square



%annotation('arrow',[0.85 0.8],[0.75 0.68])
%annotation('arrow',[0.9 0.85],[0.75 0.68])
%annotation('arrow',[0.9 0.85],[0.68 0.61])


set(gca, 'FontSize',12)

fontname(gcf, 'CMU Serif')

return;

%% 

figure(5);
s = pcolor(X,Y,H);
s.EdgeColor = 'None';
colormap(gray)
%title('H')
axis off


%% 

figure(6);
pc = pcolor(X,Y,abs(mu));
clim([0 max(abs(mu(:)))/3])
pc.EdgeColor = 'none';
%title('|\rho|')
axis off

%% 

figure(7);
pc = pcolor(X,Y,abs(real(phi_tot)));
%clim([0 max(abs(phi_n_tot(:)))*0.9])
ylim([-1500 1500])
pc.EdgeColor = 'none';
colorbar
title('Shelf displacement |Re(\phi_z)|') 
