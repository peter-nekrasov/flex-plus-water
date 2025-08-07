%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for the 
% Helmholtz scattering problem:
%
%      \Delta u + k^2 (1 + V(x) ) u = f
%
% in this case, f = - k^2 V exp(i k x)
%
%%%%%

L = 500;
N = 201; % needs to be an odd number

zk = 0.1;

xs = L*(-floor(N/2):floor(N/2))/floor(N/2);
[xxgrid,yygrid] = meshgrid(xs);

h = xs(2) - xs(1);

coefs = bump2_helm(xxgrid,yygrid,-0.5,50);
V = coefs{1};

dinds = find(abs(V) > 1e-12 );
[iinds,jinds] = find(abs(V) > 1e-12 );

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

[src,targ,ind,sz] = get_fft_grid(N,L);
[inds,corrs] = get_correct_helm(h);
spmats = get_sparse_corr(size(xxgrid),inds,corrs);
kerns = kernmat(src,targ,@(s,t) helm2d.green_cell_helm(zk,s,t),h);
kerns = gen_fft_kerns2(kerns,sz,ind);

% Solve with GMRES
start = tic;
sol = gmres(@(mu) fast_apply_fft_sub_helm(mu,kerns,zk,coefs,spmats,h,dinds,iinds,jinds,xxgrid),rhs_vec,[],1e-12,200);
mu = zeros(size(xxgrid));
mu(dinds) = sol;
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

evalkerns = {kerns{1}};
evalcorrs = {spmats{1}};

usca = sol_eval_fft_sub_helm(sol,evalkerns,evalcorrs,h,dinds,iinds,jinds,xxgrid);

utot = usca + uinc;

%%

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

%%

figure(1);
s = surf(xxgrid,yygrid,H);
s.EdgeColor = 'none';

figure(2);
s = surf(xxgrid,yygrid,real(phi_n));
s.EdgeColor = 'none';
