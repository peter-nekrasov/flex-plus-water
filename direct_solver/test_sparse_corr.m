%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for plane wave 
% scattering of flexural-gravity waves
%
%
%%%%%

clear 
close all
addpath(genpath('..'))

L = 500;
h = 5;

xs = -L:h:L;
xl = -2*L:h:2*L;
[~,n] = size(xs);
[X,Y] = meshgrid(xs);
[XL,YL] = meshgrid(xl);
w = 8;

[coefs, H] = bump2(X,Y,2,50,w); % remove gbar from coefs vector
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

% mu = zeros(size(X));
mu = exp(-X.^2/1000 - Y.^2/1000);
% mu(100:101,101) = 1;
mu = mu(:);

% Constructing full fft original integral operators
[inds,corrs] = get_correct(h,a0);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
ind = find((XL == 0) & (YL ==0));
sz = size(XL);

kerns = gen_fft_kerns(kerns,sz,ind);

[v1,Gs_mu] = fast_apply_fft(mu,kerns,coefs);

% Constructing full fft original integral operators
spmat = get_sparse_corr(size(X),inds,corrs);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
ind = find((XL == 0) & (YL ==0));
sz = size(XL);

kerns = gen_fft_kerns(kerns,sz,ind);

[v2,Gs_mu2] = fast_apply_fft_plus_corr(mu,kerns,coefs,spmat,h);

norm(Gs_mu(:) - Gs_mu2(:))

err = abs(v1-v2)

figure(1); clf;
tiledlayout(2,2,'TileSpacing','compact');
nexttile 
pcolor(X,Y,imag(reshape(v1, size(X))),'EdgeColor','none')
nexttile
pcolor(X,Y,imag(reshape(v2, size(X))),'EdgeColor','none')
nexttile
pcolor(X,Y,(reshape(err,size(X))),'EdgeColor','none')
colorbar
nexttile
pcolor(X,Y,(reshape(mu,size(X))),'EdgeColor','none')

return

figure(2);
tiledlayout(2,3)

nexttile
pc = pcolor(X,Y,real(mu));
pc.EdgeColor = 'none';
title('Re(\mu)')
colorbar

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
err = get_fin_diff_err(X,Y,mu,phi_n_tot,phi_tot,h,coefs,10,10)

return

%%

figure(1);
s = surf(X,Y,H);
s.EdgeColor = 'none';

figure(2);
s = surf(X,Y,real(phi_n));
s.EdgeColor = 'none';
