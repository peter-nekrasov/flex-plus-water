%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for the 
% Helmholtz scattering problem:
%
%      \Delta u + k^2 (1 + V(x) ) u = f
%
% in this case, f = - k^2 V exp(i k x)
%
% Solved iteratively using FFT + GMRES 
%
%%%%%

L = 5;
Ns = 1 + 50:50:700; % needs to be an odd number for FFT

zk = 4;

errs = zeros(length(Ns),1);

for ii = 1:length(Ns)

N = Ns(ii);
xs = L*(-floor(N/2):floor(N/2))/floor(N/2);
[xxgrid,yygrid] = meshgrid(xs);

h = xs(2) - xs(1);

coefs = bump2_helm(xxgrid,yygrid,-0.5,0.5);
V = coefs{1};

dinds = find(abs(V) > 1e-12 );
[iinds,jinds] = find(abs(V) > 1e-12 );

% RHS (Incident field)
k1 = zk;
k2 = 0;
uinc = exp(1i*k1*xxgrid+1i*k2*yygrid);
[rhs_vec, rhs] = get_rhs_vec_helm(coefs,zk,uinc);
rhs_vec = rhs_vec(dinds);

% Constructing integral operators

[src,targ,ind,sz] = get_fft_grid(N,L);
[inds,corrs] = get_correct_helm(h,zk);
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

errs(ii) = get_fin_diff_err_helm(xxgrid,yygrid,utot,h,coefs,0.1,0.1,zk);

end

%%

figure(2);
plot(log10(Ns),log10(errs),'x-');
hold on
plot(log10(Ns),log10(1e9*Ns.^(-8)))
xlabel('N')
legend('error','N^{-8}')