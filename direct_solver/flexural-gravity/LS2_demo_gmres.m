%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for the 
% flexural-gravity wave scattering problem
%
% Solved iteratively using FFT + GMRES
%
%%%%%

L = 500;
N = 405; % needs to be an odd number

xs = L*(-floor(N/2):floor(N/2))/floor(N/2);
[xxgrid,yygrid] = meshgrid(xs);

w = 8;

h = xs(2) - xs(1);

[coefs, H] = bump2(xxgrid,yygrid,2,50,w); % remove gbar from coefs vector
E = 7E9;

a0 = coefs{1}; 
abar = coefs{2};
b0 = coefs{3}; 
g0 = coefs{5}; 

dinds = find((abar / a0) > 1e-12 );
[iinds,jinds] = find((abar / a0) > 1e-12 );

% Finding positive real roots
[rts,ejs] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0));
ejs = ejs/a0;

% RHS (Incident field)
k1 = k;
k2 = 0;
phiinc = exp(1i*k1*xxgrid+1i*k2*yygrid);
[rhs_vec, rhs] = get_rhs_vec(coefs,k1,k2,phiinc);
rhs_vec = rhs_vec(dinds) / k;
rhs = rhs / k;

figure(1);
tiledlayout(1,4);

nexttile
s = pcolor(xxgrid,yygrid,H);
s.EdgeColor = 'None';
colorbar
title('H')
drawnow

nexttile
s = pcolor(xxgrid,yygrid,(coefs{1} + coefs{2}));
s.EdgeColor = 'None';
colorbar
title('\alpha')
drawnow

nexttile
s = pcolor(xxgrid,yygrid,(coefs{3} + coefs{4}));
s.EdgeColor = 'None';
colorbar
title('\beta')
drawnow

nexttile
s = pcolor(xxgrid,yygrid,real(rhs));
s.EdgeColor = 'None';
colorbar
title('rhs')
drawnow

% Constructing integral operators

[src,targ,ind,sz] = get_fft_grid(N,L);
[inds,corrs] = get_correct(h,a0);
spmat = get_sparse_corr(size(xxgrid),inds,corrs);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
kerns = gen_fft_kerns(kerns,sz,ind);

% Solve with GMRES
start = tic;
sol = gmres(@(mu) fast_apply_fft_sub(mu,kerns,coefs,spmat,h,dinds,iinds,jinds,xxgrid),rhs_vec,[],1e-12,200);
mu = zeros(size(xxgrid));
mu(dinds) = sol;
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

evalkerns = {kerns{1}, kerns{4}};
evalcorrs = {spmat{1}, spmat{4}};

[phi, phi_n] = sol_eval_fft_sub(sol,evalkerns,evalcorrs,h,dinds,iinds,jinds,xxgrid);

phi_tot = phi + phiinc/k;
phi_n_tot = phi_n + phiinc;

%%

figure(2);
tiledlayout(2,3)

nexttile
pc = pcolor(xxgrid,yygrid,real(mu));
pc.EdgeColor = 'none';
title('Re(\mu)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,real(phi_tot));
pc.EdgeColor = 'none';
title('Re(\phi)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(phi_tot));
pc.EdgeColor = 'none';
title('|\phi|')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,real(phi_n_tot));
pc.EdgeColor = 'none';
title('real(\phi_n)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(phi_n_tot));
pc.EdgeColor = 'none';
title('|\phi_n|')
colorbar
       
% Calculate error with finite difference
err = get_fin_diff_err(xxgrid,yygrid,mu,phi_n_tot,phi_tot,h,coefs,10,10)

return

%%

figure(2);
tiledlayout(1,2)

nexttile
pc = pcolor(xxgrid,yygrid,real(phi_n));
pc.EdgeColor = 'none';
title('real(\phi_n)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,real(phi));
pc.EdgeColor = 'none';
title('|\phi_n|')
colorbar