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
h = 50;

xs = -L:h:L;
xl = -2*L:h:2*L;
[~,n] = size(xs);
[X,Y] = meshgrid(xs);
[XL,YL] = meshgrid(xl);

[coefs, H] = bumps(X,Y,-8000,8000,0.4,100); % remove gbar from coefs vector
E = 7E9;

a0 = coefs{1}; 
b0 = coefs{3}; 
g0 = coefs{5}; 

% Finding positive real roots
[rts,ejs] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0));
ejs = ejs/a0;

src = [0;0];
targ = [XL(:).'; YL(:).'];

% RHS (Incident field)
k1 = k*cos(pi/4);
k2 = k*sin(pi/4);
phiinc = exp(1i*k1*X+1i*k2*Y);
[rhs_vec, rhs] = get_rhs_vec(coefs,k1,k2,phiinc);

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
s = pcolor(X,Y,real(E*rhs));
s.EdgeColor = 'None';
colorbar
title('rhs')
drawnow


% Constructing integral operators
[inds,corrs] = get_correct(h,a0);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs), h,inds,corrs);

ind = find((XL == 0) & (YL ==0));
sz = size(XL);

kerns = gen_fft_kerns(kerns,sz,ind);

evalkerns = {kerns{1}, kerns{4}};

% Solve with GMRES
start = tic;
mu2 = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,[],1e-12,200);
mu = reshape(mu2, size(X));
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

[phi, phi_n] = sol_eval_fft(mu,evalkerns);

phi_tot = phi + phiinc;
phi_n_tot = phi_n + k*phiinc;

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
err = get_fin_diff_err(X,Y,mu,phi_n_tot,phi_tot,h,coefs,200,200)

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
pc = pcolor(X,Y,real(phi_n_tot));
pc.EdgeColor = 'none';
title('Real(\phi_n)')
colorbar
clim([-0.06 0.06])

%%
V = coefs;
n = numel(V{2});
slf = 1:n;
a0 = V{1};
abar = V{2}(slf);
b0 = V{3};
bbar = V{4}(slf);
g0 = V{5};
gbar = V{6};
alphax = V{7}(slf);
alphay = V{8}(slf);
alphaxx = V{9}(slf);
alphaxy = V{10}(slf);
alphayy = V{11}(slf);
nu = V{end};


[inds, corrs] = get_correct(h,a0);
scorr = get_sparse_corr(size(X),inds,corrs);

rs = [X(:).';Y(:).'];

n = size(rs,2);

Afun = @(i,j) kernbyindex(i,j,rs,h,scorr,rts,ejs,coefs);
rank_or_tol = 1E-8;
pxyfun = [];
opts = [];
opts.verb = true;
occ = 500;

%%
rprop = rts(abs(imag(rts))<1E-12);
nwave = (1.5*L/4)*rprop;
width = 1;
q     = 16;
p     = ceil(10*nwave);
p0    = 50;
[proxy,pw] = proxy_ann_pts2(p,p0,q,width);

t1 = max(abs((a0*bbar-abar*b0)./a0));
t2 = max(abs((0.5*abar.*g0./a0)));

pxfun = @(x,slf,nbr,l,ctr) proxyfun(slf,nbr,l,ctr,x,h,rts,ejs,V,proxy,pw,t1,t2);

%%
Afun = @(i,j) kernbyindex(i,j,rs,h,scorr,rts,ejs,coefs);
rank_or_tol = 1E-6;
opts = [];
opts.verb = true;
occ = 2000;
F2 = rskelf(Afun,rs,occ,rank_or_tol,pxfun,opts)

%%
sol2 = rskelf_sv(F2,rhs_vec);
