%%%%%
%
% Continuous scattering with variable beta 
%
%%%%%

addpath('../greens function/')

% Parameters
h = 50;
xs = -5000:h:5000;
[~,n] = size(xs);
[X,Y] = meshgrid(xs);
H0 = 20;
w = 1;
a0 = 6.41026E8*H0^3;
b0 = (917*H0*w^2 - 9800) / a0;
g0 = - 1000*w^2 / a0;

% Finding positive real roots
[rts,~] = find_roots(b0,g0);
k = rts((imag(rts) == 0) & (real(rts) > 0));

% Perturbing coefficients
geo = gaussian(X,Y,H0,5/k);
H = geo{1};
beta = (917*H*w^2 - 9800) ./ geo{2};
gamma = -1000*w^2 ./ geo{2};

% Constructing integral operator
gf = green(X - min(xs), Y - min(xs), beta0, gamma, false);
Gs = gf{1};
% fill in diagonal correction 
Gs_aug = [Gs, flip(Gs(1:end,2:end),2); ...
    flip(Gs(2:end,1:end)), flip(flip(Gs(2:end,2:end)),2)];
Gs_aug_hat = fft2(Gs_aug)*h*h;

% RHS (Incident field)
phiinc = exp(1i*k*X);
rhs = k*betabar.*phiinc;
rhs_vec = rhs(:);

% Solve with GMRES
start = tic;
mu = gmres(@(mu) lhs(mu,betabar,Gs_aug_hat),rhs_vec,[],1e-12,200);
mu = reshape(mu, size(X));
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

% Evaluation and plotting
gf = green(X - min(xs), Y - min(xs), beta0, gamma, true);
Gc = gf{1};
% fill in diagonal correction 
Gc_aug = [Gc, flip(Gc(1:end,2:end),2); ...
    flip(Gc(2:end,1:end)), flip(flip(Gc(2:end,2:end)),2)];
Gc_aug_hat = fft2(Gc_aug)*h*h;
mu_aug = [mu, zeros(n,n-1); zeros(n-1,n), zeros(n-1)];
mu_aug_hat = fft2(mu_aug);
phi_aug = ifft2(Gc_aug_hat.*mu_aug_hat);
phi_n_aug = ifft2(Gs_aug_hat.*mu_aug_hat);
phi = phi_aug(1:n,1:n);
phi_n = phi_n_aug(1:n,1:n);

phi_tot = phiinc+phi;
phi_n_tot = phi_n + k*phiinc;

tiledlayout(1,5)
nexttile
pc = pcolor(X,Y,beta);
pc.EdgeColor = 'none';
title('\beta')
colorbar

nexttile
pc = pcolor(X,Y,real(phi_tot));
pc.EdgeColor = 'none';
%clim([-1.5 1.5])
title('Re(\phi)')
colorbar

nexttile
pc = pcolor(X,Y,abs(phi_tot));
pc.EdgeColor = 'none';
%clim([0 1.5])
title('|\phi|')
colorbar

nexttile
pc = pcolor(X,Y,real(phi_n_tot));
pc.EdgeColor = 'none';
%clim([-1.5 1.5]*k)
title('real(\phi_n)')
colorbar

nexttile
pc = pcolor(X,Y,abs(phi_n_tot));
pc.EdgeColor = 'none';
%clim([0 1.5]*k)
title('|\phi_n|')
colorbar

rmpath('../greens function/')

%% Calculating error

ind = intersect(find(X == 16), find(Y == 0));
[ii, jj] = ind2sub(size(X),ind);
% disp(phi(ii,jj))

% finite difference stencils
d1 = zeros(9,1);
d1(1) = 7/240;	
d1(2) = -2/5;
d1(3) = 169/60;
d1(4) = -122/15;
d1(5) = 91/8;
d1(6) = -122/15;
d1(7) = 169/60;
d1(8) = -2/5;
d1(9) = 7/240;

d2 = zeros(9, 1);
d2(1) = -1/560;
d2(2) = 8/315;
d2(3) = -1/5;
d2(4) = 8/5;
d2(5) = -205/72;
d2(6) = 8/5;
d2(7) = -1/5;
d2(8) = 8/315;
d2(9) = -1/560;

bilap = zeros(9,9);
bilap(:,5) = d1;
bilap(5,:) = bilap(5,:) + d1.';
bilap = bilap + 2*(d2*d2.');
bilap = bilap / h^4;

phi_n_tot_sub = phi_n_tot(ii-4:ii+4,jj-4:jj+4);
first = sum(bilap.*phi_n_tot_sub,'all') ;
second = -beta(ii,jj)*phi_n_tot(ii,jj);
third = gamma*phi_tot(ii,jj);
err = abs(first + second + third) / (abs(phi_n(ii,jj)))

% phiinc_sub = phiinc(ii-4:ii+4,jj-4:jj+4);
% first = k*sum(bilap.*phiinc_sub,'all') / h^4;
% second = -beta0*k*phiinc(ii,jj);
% third = gamma*phiinc(ii,jj);
% err = first + second + third

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