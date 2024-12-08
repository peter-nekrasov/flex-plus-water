%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for plane wave 
% scattering of flexural-gravity waves with random gaussian thickness
%
%%%%%

clear 
close all
addpath(genpath('..'))

h = 20;

x1 = -2E3;
x2 = 6E3;

y1 = -4E3;
y2 = 4E3;

xs = x1:h:x2;
ys = y1:h:y2;
xl = 2*x1:h:2*x2;
yl = 2*y1:h:2*y2;
[~,n] = size(xs);
[X,Y] = meshgrid(xs,ys);
[XL,YL] = meshgrid(xl,yl);

[coefs, H] = rolls(X,Y,0,50E2,-30E2,30E2,1,500,0.5); % remove gbar from coefs vector

a0 = coefs{1}; 
b0 = coefs{3}; 
g0 = coefs{5}; 

% Finding positive real roots
[rts,ejs] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0));
ejs = ejs/a0;

% RHS (Incident field)
% k1 = k*cos(0);
% k2 = k*sin(0);
% phiinc = exp(1i*k1*X+1i*k2*Y);
% phininc = k*exp(1i*k1*X+1i*k2*Y);
% [rhs_vec, rhsp] = get_rhs_vec(coefs,k1,k2,phiinc);

src = [-x1/2;0]+1i*[3000;0];
targ = [X(:).'; Y(:).'];

% RHS (Incident field)
kerns = green(src,targ,rts,ejs);
cnst = max(kerns{1}(:));
kerns = cellfun(@(x) x/cnst,kerns,'UniformOutput',false);
phininc = reshape(kerns{1},size(X));
phiinc = reshape(kerns{4},size(X));
[rhs_vec] = get_rhs_vec2(coefs,kerns);
rhsp = reshape(rhs_vec,size(X));

figure(1);
tiledlayout(1,5);

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

nexttile
s = pcolor(X,Y,real(phininc));
s.EdgeColor = 'None';
colorbar
title('phi_n incident')
drawnow



% Constructing integral operators
src = [xl(ceil(end/2)); yl(ceil(end/2))];
targ = [XL(:).'; YL(:).'];

[inds,corrs] = get_correct(h,a0);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs), h, inds,corrs);

ind = find((XL == src(1)) & (YL == src(2)));
sz = size(XL);

kerns = gen_fft_kerns(kerns,sz,ind);

evalkerns = {kerns{1}, kerns{4}};

% Solve with GMRES
start = tic;
mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,20,1e-9,1000);
mu = reshape(mu, size(X));
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

[phi, phi_n] = sol_eval_fft(mu,evalkerns);

phi_tot = phi + phiinc;
phi_n_tot = phi_n + phininc;

%%


figure(2);
pc = pcolor(X,Y,real(mu));
pc.EdgeColor = 'none';
title('Real(\mu)')
colorbar

figure(3);
tiledlayout(2,2)

nexttile
pc = pcolor(X,Y,real(phi_tot));
pc.EdgeColor = 'none';
title('Re(\phi)')
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
colorbar

nexttile
pc = pcolor(X,Y,abs(phi_n_tot));
pc.EdgeColor = 'none';
title('|\phi_n|')
colorbar
       
% Calculate error with finite difference
err = get_fin_diff_err(X,Y,mu,phi_n_tot,phi_tot,h,coefs,1000,400)


return

%% Figure generation for Jeremy

figure(4);

t = tiledlayout(1,4,"TileSpacing","tight");

X1 = X / 1000 + 2;
Y1 = Y / 1000 + 5;

ax1 = nexttile
s = pcolor(X1,Y1,H);
s.EdgeColor = 'None';
colormap(ax1,gray)
clim([min(H(:)) max(H(:))])
colorbar
title('H')


nexttile
pc = pcolor(X1,Y1,abs(mu));
clim([0 max(abs(mu(:)))])
pc.EdgeColor = 'none';
colorbar
title('|\rho|')

nexttile
pc = pcolor(X1,Y1,real(phi_tot));
pc.EdgeColor = 'none';
title('Re(\phi)')
colorbar

nexttile
pc = pcolor(X1,Y1,abs(phi_tot));
pc.EdgeColor = 'none';
title('|\phi|')
colorbar

%% Figure generation for Jeremy

figure(4);

t = tiledlayout('flow','TileSpacing','tight'); 

X1 = X / 1000 + 5;
Y1 = Y / 1000 + 5;

ax1 = nexttile
s = pcolor(X1,Y1,H*0+3);
s.EdgeColor = 'None';
colormap(ax1,gray)
clim([min(H(:)) max(H(:))])
colorbar
title('H')


nexttile
pc = pcolor(X1,Y1,H*0);
clim([0 max(abs(mu(:)))/3])
pc.EdgeColor = 'none';
colorbar
title('|\rho|')

nexttile([2 2]);
pc = pcolor(X1,Y1,H*0);
clim([0 max(abs(phi_n_tot(:)))*0.75])
pc.EdgeColor = 'none';
colorbar
title('|\phi_n|')


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
pc = pcolor(X,Y,abs(phi_n_tot));
clim([0 max(abs(phi_n_tot(:)))*0.75])
pc.EdgeColor = 'none';
%title('|\phi_n|')
axis off

