h = 0.5;
L = 70;

xs = -L/2:h:L/2;
xl = -L:h:L;

ys = -L/2:h:L/2;
yl = -L:h:L;

[~,~,coefs,Hs] = get_Y(-35,35,-5,65,h,-0.98,200);

[X,Y] = meshgrid(xs,ys);
[XL,YL] = meshgrid(xl,yl);

a0 = coefs{1}; 
b0 = coefs{3}; 
g0 = coefs{5}; 

% Finding positive real roots
[rts,ejs] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0));
ejs = ejs/a0;

[inds, corrs] = get_correct(h,a0);
scorr = get_sparse_corr(size(X),inds,corrs);

rs = [X(:).';Y(:).'];

%%
src = [0;0];
targ = [XL(:).'; YL(:).'];

% RHS (Incident field)
k1 = k*cos(pi/3);
k2 = k*sin(pi/3);
phiinc = exp(1i*k1*X+1i*k2*Y);
[rhs_vec, rhs] = get_rhs_vec(coefs,k1,k2,phiinc);

% Constructing integral operators
[inds,corrs] = get_correct(h,a0);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h, inds,corrs);

ind = find((XL == 0) & (YL ==0));
sz = size(XL);

kerns = gen_fft_kerns(kerns,sz,ind);

evalkerns = {kerns{1}, kerns{4}};

%%
rprop = rts(abs(imag(rts))<1E-12);
nwave = (1.5*L/4)*rprop;
width = 1;
q     = 16;
p     = ceil(10*nwave);
[proxy,pw] = proxy_ann_pts(p,q,width);

abar = coefs{2};
bbar = coefs{4};

t1 = max(abs((a0*bbar(:)-abar(:)*b0)./a0));
t2 = max(abs((0.5*abar(:).*g0./a0)));

pxfun = @(x,slf,nbr,l,ctr) proxyfun(slf,nbr,l,ctr,x,h,rts,ejs,coefs,proxy,pw,t1,t2);

%%
Afun = @(i,j) kernbyindex(i,j,rs,h,scorr,rts,ejs,coefs);
rank_or_tol = 1E-1;
opts = [];
opts.verb = true;
occ = 500;
F2 = rskelf(Afun,rs,occ,rank_or_tol,pxfun,opts)

%%
sol2 = rskelf_sv(F2,rhs_vec);

mu = gmres(@(mu) fast_apply_fft_w_preconditioner(mu,kerns,coefs,F2),sol2,150,1e-8,1000);


%%

%mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,100,1e-9,200);

load gong.mat
sound(y)


%%

mu = reshape(mu,size(X));
sol2 = reshape(sol2,size(X));

figure(1);
tiledlayout(1,2)
nexttile
pl1 = pcolor(X,Y,abs(mu));
pl1.EdgeColor = 'none';
colorbar

nexttile
pl2 = pcolor(X,Y,abs(sol2));
pl2.EdgeColor = 'none';
colorbar

%%

[phi, phi_n] = sol_eval_fft(sol2,evalkerns);

phiinc = exp(1i*k1*X+1i*k2*Y);
phininc = k*exp(1i*k1*X+1i*k2*Y);

phi_tot = phi + phiinc;
phi_n_tot = phi_n + phininc;

figure(2);
tiledlayout(1,2)
nexttile
pl1 = pcolor(X,Y,real(phi_tot.'));
pl1.EdgeColor = 'none';
colorbar

nexttile
pl2 = pcolor(X,Y,real(phi_n_tot.'));
pl2.EdgeColor = 'none';
colorbar