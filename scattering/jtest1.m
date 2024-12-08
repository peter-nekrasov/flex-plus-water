L = 2000;
h = 50;

xs = -L:h:L;
xl = -2*L:h:2*L;
[~,n] = size(xs);
[X,Y] = meshgrid(xs);
[XL,YL] = meshgrid(xl);

[coefs, H] = bump2(X,Y,5,200); % remove gbar from coefs vector
E = 7E9;

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
k1 = 2*sqrt(2)*k/3;
k2 = k/3;
phiinc = exp(1i*k1*X+1i*k2*Y);
[rhs_vec, rhs] = get_rhs_vec(coefs,k1,k2,phiinc);

% Constructing integral operators
[inds,corrs] = get_correct(h,a0);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h, inds,corrs);

ind = find((XL == 0) & (YL ==0));
sz = size(XL);

kerns = gen_fft_kerns(kerns,sz,ind);

evalkerns = {kerns{1}, kerns{4}};

rhs_vec = zeros(size(rs,2),1);
ind = 3281;
rhs_vec(ind) = 1;
vout = fast_apply_fft(rhs_vec,kerns,coefs);
% Solve with GMRES
start = tic;
%mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,[],1e-12,200);

%%
n = size(rs,2);
[xmat] = kernbyindex(1:n,ind,rs,h,scorr,rts,ejs,coefs);

%%
Afun = @(i,j) kernbyindex(i,j,rs,h,scorr,rts,ejs,coefs);
rank_or_tol = 1E-8;
pxyfun = [];
opts = [];
opts.verb = true;
occ = 500;
F = rskelf(Afun,rs,occ,rank_or_tol,pxyfun,opts)

%%
% RHS (Incident field)
k1 = 2*sqrt(2)*k/3;
k2 = k/3;
phiinc = exp(1i*k1*X+1i*k2*Y);
[rhs_vec, rhs] = get_rhs_vec(coefs,k1,k2,phiinc);

%%
sol = rskelf_sv(F,rhs_vec);

%%
mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,[],1e-12,200);

%%
rprop = rts(abs(imag(rts))<1E-12);
nwave = (1.5*L/4)*rprop;
width = 1;
q     = 16;
p     = ceil(10*nwave);
[proxy,pw] = proxy_ann_pts(p,q,width);

t1 = max(abs((a0*bbar-abar*b0)./a0));
t2 = max(abs((0.5*abar.*g0./a0)));

pxfun = @(x,slf,nbr,l,ctr) proxyfun(slf,nbr,l,ctr,x,h,rts,ejs,V,proxy,pw,t1,t2);

%%
Afun = @(i,j) kernbyindex(i,j,rs,h,scorr,rts,ejs,coefs);
rank_or_tol = 1E-8;
opts = [];
opts.verb = true;
occ = 500;
F2 = rskelf(Afun,rs,occ,rank_or_tol,pxfun,opts)

%%
sol2 = rskelf_sv(F2,rhs_vec);

%%
mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,[],1e-12,200);