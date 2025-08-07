%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for the 
% Helmholtz scattering problem:
%
%      \Delta u + k^2 (1 + V(x) ) u = f
%
% in this case, f = - k^2 V exp(i k x)
%
% Solved directly using skeletonization
%
%%%%%

L = 500;
N = 201; 

zk = 0.1;

xs = L*(-floor(N/2):floor(N/2))/floor(N/2);
[xxgrid,yygrid] = meshgrid(xs);

h = xs(2) - xs(1);

coefs = bump2_helm(xxgrid,yygrid,-0.5,50);
V = coefs{1};

dinds = find(abs(V) > 1e-12 );
[iinds,jinds] = find(abs(V) > 1e-12 );

srcinfo = []; srcinfo.r = [xxgrid(dinds) yygrid(dinds)].'; srcinfo.wts = h^2*ones(length(dinds),1);
targinfo = []; targinfo.r = [xxgrid(dinds) yygrid(dinds)].'; 
targinfo.V = V(dinds);

% RHS (Incident field)
k1 = zk;
k2 = 0;
uinc = exp(1i*k1*xxgrid+1i*k2*yygrid);
[rhs_vec, rhs] = get_rhs_vec_helm(coefs,zk,uinc);
rhs_vec = rhs_vec(dinds);

figure(1); clf
tiledlayout(1,2);

nexttile
s = pcolor(xxgrid,yygrid,V);
s.EdgeColor = 'None';
colorbar
title('V')
drawnow

nexttile
s = pcolor(xxgrid,yygrid,real(rhs));
s.EdgeColor = 'None';
colorbar
title('rhs')
drawnow

% Constructing integral operators

[inds,corrs] = get_correct_helm(h);
spmats = get_sparse_corr(size(xxgrid),inds,corrs);

% Constructing integral operators

[inds,corrs] = get_correct_helm(h);
spmat = get_sparse_corr(size(xxgrid),inds,corrs);
idspmat = id_plus_corr_sum_helm(zk,coefs,spmat,dinds,h);
kernfun = @(s,t) kern_sum_helm(zk,s,t);
Afun = @(i,j) kern_matgen(i,j,srcinfo,targinfo,idspmat,kernfun);

% Solve with FLAM

x = srcinfo.r;
occ = 2000;
rank_or_tol = 1e-8;
pxyfun = [];
opts = [];

start = tic;
F = rskelf(Afun,x,occ,rank_or_tol,pxyfun,opts);
t2 = toc(start);
fprintf('%5.2e s : time to factorize inverse \n',t2)

start = tic;
sol2 = rskelf_sv(F,rhs_vec);
t3 = toc(start);
fprintf('%5.2e s : time to solve (skel) \n',t3)

mu = zeros(size(xxgrid));
mu(dinds) = sol2;

% Plot with FFT

[src,targ,ind,sz] = get_fft_grid(N,L);
kerns = kernmat(src,targ,@(s,t) helm2d.green_cell_helm(zk,s,t),h);
kerns = gen_fft_kerns2(kerns,sz,ind);

evalkerns = {kerns{1}};
evalcorrs = {spmat{1}};

usca = sol_eval_fft_sub_helm(sol2,evalkerns,evalcorrs,h,dinds,iinds,jinds,xxgrid);

utot = usca + uinc;

figure(2);
tiledlayout(1,3)

nexttile
pc = pcolor(xxgrid,yygrid,real(mu));
pc.EdgeColor = 'none';
title('Re(\mu)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,real(utot));
pc.EdgeColor = 'none';
title('Re(\phi)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(utot));
pc.EdgeColor = 'none';
title('|\phi|')
colorbar
       
% Calculate error with finite difference
err = get_fin_diff_err_helm(xxgrid,yygrid,utot,h,coefs,10,10,zk)

return
