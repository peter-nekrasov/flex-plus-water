addpath('greens function/')

% Parameters
h = 0.5;
xs = -50:h:50;
[~,n] = size(xs);
[X,Y] = meshgrid(xs);
beta0 = 5;
gamma = -1;

% Finding positive real roots
[rts,~] = find_roots(beta0,gamma);
k = rts((imag(rts) == 0) & (real(rts) > 0));

% Perturbing coefficients
betabar = 1*exp(-(X.^2 + Y.^2)/25); 
beta = beta0 + betabar;

% Constructing integral operator
Gs = green(sqrt((X - min(xs)).^2 + (Y- min(xs)).^2), beta0, gamma, false);
Gs(1,1) = 0; % fill in diagonal correction 
Gs_aug = [Gs, flip(Gs(1:end,2:end),2); ...
    flip(Gs(2:end,1:end)), flip(flip(Gs(2:end,2:end)),2)];
Gs_aug_hat = fft2(Gs_aug)*h*h;


% RHS (Incident field)
phiinc = exp(1i*k*X);
rhs = k*betabar.*phiinc;
rhs_vec = rhs(:);

% Solve with GMRES
start = tic;
mu = gmres(@(mu) lhs(mu,betabar,Gs_aug_hat),rhs_vec,[],1e-13,200);
mu = reshape(mu, size(X));
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

% Evaluation and plotting
Gc = green(sqrt((X - min(xs)).^2 + (Y- min(xs)).^2), beta0, gamma, true);
Gc(1,1) = 0; % fill in diagonal correction 
Gc_aug = [Gc, flip(Gc(1:end,2:end),2); ...
    flip(Gc(2:end,1:end)), flip(flip(Gc(2:end,2:end)),2)];
Gc_aug_hat = fft2(Gc_aug)*h*h;
mu_aug = [mu, zeros(n,n-1); zeros(n-1,n), zeros(n-1)];
mu_aug_hat = fft2(mu_aug);
phi_aug = ifft2(Gc_aug_hat.*mu_aug_hat);
phi_n_aug = ifft2(Gs_aug_hat.*mu_aug_hat);
phi = phi_aug(1:n,1:n);
phi_n = phi_n_aug(1:n,1:n);

tiledlayout(1,5)
nexttile
h = pcolor(X,Y,beta);
h.EdgeColor = 'none';
title('\beta')
colorbar

nexttile
h = pcolor(X,Y,real(phiinc+phi));
h.EdgeColor = 'none';
clim([-1.5 1.5])
title('Re(\phi)')
colorbar

nexttile
h = pcolor(X,Y,abs(phiinc+phi));
h.EdgeColor = 'none';
clim([0 1.5])
title('|\phi|')
colorbar

nexttile
h = pcolor(X,Y,real(phi_n + k*phiinc));
h.EdgeColor = 'none';
clim([-1.5 1.5]*k)
title('real(\phi_n)')
colorbar

nexttile
h = pcolor(X,Y,abs(phi_n + k*phiinc));
h.EdgeColor = 'none';
clim([0 1.5]*k)
title('|\phi_n|')
colorbar


ind = intersect(find(X == 16), find(Y == 0));
disp(phi(ind))

rmpath('greens function/')

function v = lhs(mu,betabar,Gs_aug_hat)
    N = sqrt(size(mu));
    N = N(1);
    mu = reshape(mu, [N N]);
    mu_aug = [mu, zeros(N,N-1); zeros(N-1,N), zeros(N-1)];
    mu_aug_hat = fft2(mu_aug);
    Gs_mu_aug = ifft2(Gs_aug_hat.*mu_aug_hat);
    Gs_mu = Gs_mu_aug(1:N, 1:N);
    v = mu - betabar.*Gs_mu;
    v = v(:);
end