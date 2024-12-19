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

L = 750;
h = 1;
H0 = 1;

xs = -L:h:L-h;
xl = -2*L:h:2*L-2*h;
[X,Y] = meshgrid(xs);
[XL,YL] = meshgrid(xl);

% xeval = -1.5*L:h:2*L;
% xleval = -3*L:h:4*L;
% [Xeval,Yeval] = meshgrid(xeval);
% [XLeval,YLeval] = meshgrid(xleval);

load('cracks.mat');
nc = -.995;
foutx = nc*(foutx(1:h:end,1:h:end))/max(fout(:));
fouty = nc*(fouty(1:h:end,1:h:end))/max(fout(:));
foutxx = nc*(foutxx(1:h:end,1:h:end))/max(fout(:));
foutxy = nc*(foutxy(1:h:end,1:h:end))/max(fout(:));
foutyy = nc*(foutyy(1:h:end,1:h:end))/max(fout(:));
fout = nc*(fout(1:h:end,1:h:end)-min(fout(:)))/max(fout(:))+H0;



figure(1);
tiledlayout(2,3);

nexttile
s = pcolor(X,Y,fout);
s.EdgeColor = 'None';
colorbar
title('H')
drawnow

nexttile
s = pcolor(X,Y,foutx);
s.EdgeColor = 'None';
colorbar
title('\alpha')
drawnow

nexttile
s = pcolor(X,Y,fouty);
s.EdgeColor = 'None';
colorbar
title('\beta')
drawnow


nexttile
s = pcolor(X,Y,foutxx);
s.EdgeColor = 'None';
colorbar
title('rhs')
drawnow

nexttile
s = pcolor(X,Y,foutxy);
s.EdgeColor = 'None';
colorbar
title('rhs')
drawnow

nexttile
s = pcolor(X,Y,foutyy);
s.EdgeColor = 'None';
colorbar
title('rhs')
drawnow


[coefs, H] = get_coefs_from_height(fout,foutx,fouty,foutxx,foutxy,foutyy,2); % remove gbar from coefs vector
E = 7E9;

a0 = coefs{1}; 
b0 = coefs{3}; 
g0 = coefs{5}; 

% Finding positive real roots
[rts,ejs] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0))
ejs = ejs/a0;


% RHS (Incident field)
k1 = k*cos(-3*pi/4);
k2 = k*sin(-3*pi/4);
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
s = pcolor(X,Y,E*(coefs{1} + coefs{2}));
s.EdgeColor = 'None';
colorbar
title('\alpha')
drawnow

nexttile
s = pcolor(X,Y,E*(coefs{2} + coefs{3}));
s.EdgeColor = 'None';
colorbar
title('\beta')
drawnow


nexttile
s = pcolor(X,Y,real(E*rhsp));
s.EdgeColor = 'None';
colorbar
title('rhs')
drawnow


% Constructing integral operators
src = [xl(ceil(end/2)); xl(ceil(end/2))];
targ = [XL(:).'; YL(:).'];

[inds,corrs] = get_correct(h,a0);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs), h, inds,corrs);
ind = find((XL == src(1)) & (YL == src(2)));
sz = size(XL);
kerns = gen_fft_kerns(kerns,sz,ind);
ekerns = {kerns{1}, kerns{4}};

% src = [xleval(ceil(end/2)); xleval(ceil(end/2))];
% targ = [XLeval(:).'; YLeval(:).'];
% ekerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs), h, inds,corrs);
% ind = find((XLeval == src(1)) & (YLeval == src(2)));
% sz = size(XLeval);
% ekerns = gen_fft_kerns(ekerns,sz,ind);
% ekerns = {ekerns{1}, ekerns{4}};

% Solve with GMRES
start = tic;
mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,20,1e-14,300);
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

%%


figure(2);
pc = pcolor(X,Y,abs(mu));
pc.EdgeColor = 'none';
title('Abs(\mu)')
colorbar

figure(3);
tiledlayout(2,2)

nexttile
pc = pcolor(X,Y,real(phi_tot));
pc.EdgeColor = 'none';
title('Re(\phi)')
clim([-2.5 2.5])
colorbar

nexttile
pc = pcolor(X,Y,abs(phi_tot));
pc.EdgeColor = 'none';
title('|\phi|')
colorbar


nexttile
pc = pcolor(X,Y,real(phi_n_tot));
pc.EdgeColor = 'none';
title('real(\phi_n)')
%clim([-2 2])
colorbar

nexttile
pc = pcolor(X,Y,abs(phi_n_tot));
pc.EdgeColor = 'none';
title('|\phi_n|')
colorbar
       
% Calculate error with finite difference
err = get_fin_diff_err(X,Y,mu,phi_n_tot,phi_tot,h,coefs,420,-240)


return

%% Figure generation for Jeremy

f = figure(4);
f.Units = 'points';
f.InnerPosition = [600 500 600 400];

t = tiledlayout(1,2,'TileSpacing','tight'); 

X1 = X / 1000 + 5;
Y1 = Y / 1000 + 5;


c = [0.95:-0.01:0.3 ; 0.95:-0.01:0.3; 0.95:-0.01:0.3 ].';


% ax1 = nexttile;
% s = pcolor(X1,Y1,H);
% s.EdgeColor = 'None';
% colormap(ax1,c)
% clim([min(H(:)) max(H(:))])
% colorbar
% title('Thickness (m)','FontWeight','normal')


nexttile
pc = pcolor((X+L)/1000,(Y+L)/1000,abs(mu));
clim([0 5*max(abs(mu(:)))/6])
pc.EdgeColor = 'none';
xlabel('x (km)')
xlim([0 1.5])
ylim([0 1.5])
colorbar
title('|\mu|')
ylabel('y (km)')
set(gca, 'FontSize',12)
axis square




nexttile
pc = pcolor((X+L)/1000,(Y+L)/1000,abs((phi_tot)));
%clim([0 0.8*max(abs(phi_n_tot(:)))])
pc.EdgeColor = 'none';
xlim([0 1.5])
ylim([0 1.5])
colorbar
title('|\phi|')
axis square

%annotation('arrow',[0.85 0.8],[0.75 0.68])
annotation('arrow',[0.9 0.85],[0.75 0.68])
%annotation('arrow',[0.9 0.85],[0.68 0.61])


set(gca, 'FontSize',12)
xlabel('x (km)')

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
