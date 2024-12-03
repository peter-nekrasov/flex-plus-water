%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for plane wave 
% scattering of flexural-gravity waves with random gaussian thickness
%
%
%%%%%

clear 
close all
addpath(genpath('..'))

L = 10000;
h = 20;

xs = -L:h:L;
xl = -2*L:h:2*L;
[~,n] = size(xs);
[X,Y] = meshgrid(xs);
[XL,YL] = meshgrid(xl);

[coefs, H] = bumps(X,Y,-L,L,2,200); % remove gbar from coefs vector
E = 7E9;

a0 = coefs{1}; 
b0 = coefs{3}; 
g0 = coefs{5}; 

% Finding positive real roots
[rts,ejs] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0));
ejs = ejs/a0;

src = [-4000;-4000]+0*[1i;0];
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
s = pcolor(X,Y,real(phininc));
s.EdgeColor = 'None';
colorbar
title('Real(\phi^{inc}_{n})')
drawnow

nexttile
s = pcolor(X,Y,real(E*rhsp));
s.EdgeColor = 'None';
colorbar
title('rhs')
drawnow



% Constructing integral operators
src = [0; 0];
targ = [XL(:).'; YL(:).'];

[inds,corrs] = get_correct(h,a0);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs), inds,corrs);

ind = find((XL == 0) & (YL ==0));
sz = size(XL);

kerns = gen_fft_kerns(kerns,sz,ind);

evalkerns = {kerns{1}, kerns{4}};

% Solve with GMRES
start = tic;
mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,[],1e-12,200);
mu = reshape(mu, size(X));
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

[phi, phi_n] = sol_eval_fft(mu,evalkerns);

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
err = get_fin_diff_err(X,Y,mu,phi_n_tot,phi_tot,h,coefs)


%% Figure generation for Jeremy

figure(4)
s = pcolor(X,Y,H);
s.EdgeColor = 'None';
colorbar
title('H')
drawnow


figure(5) 
pc = pcolor(X,Y,abs(phi_n_tot));
pc.EdgeColor = 'none';
title('|\phi_n|')
colorbar


figure(6)
pc = pcolor(X,Y,abs(mu));
pc.EdgeColor = 'none';
title('|\rho|')
colorbar

figure(7)
pc = pcolor(X,Y,abs(real(phi_n_tot)));
pc.EdgeColor = 'none';
title('|Re(\phi_n)|')
colorbar

