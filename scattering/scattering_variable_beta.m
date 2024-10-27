%%%%%
%
% Continuous scattering with variable beta 
%
%%%%%

addpath('../greens function/')
addpath('geometry/')

% Parameters
h = 0.125;
xs = -50:h:50; %-5000:h:5000;
[~,n] = size(xs);
[X,Y] = meshgrid(xs);
% H0 = 20;
% w = 2;
a0 = 1; %6.410256410256411e+08*H0^3;
b0 = 3; %(917*H0*w^2 - 9800) ;
g0 = -1; %- 1000*w^2 ;

% Finding positive real roots
[rts,~] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0));

bbar = -X.*exp(-(X.^2 + Y.^2)/(2*(4*k)^2));
beta = b0 + bbar;


% Perturbing coefficients
% geo = gaussian(X,Y,5,H0,3*pi/k);
% H = geo{1};
% alpha = geo{2};
% beta = (917*H*w^2 - 9800) ;
% gamma = -1000*w^2 ;
% abar = alpha - a0;
% bbar = beta - b0;

% RHS (Incident field)
phiinc = exp(1i*k*X);
rhs = k*bbar.*phiinc;
rhs_vec = rhs(:);


figure(1);
tiledlayout(1,2);

% nexttile
% s = pcolor(X,Y,H);
% s.EdgeColor = 'None';
% colorbar
% title('H')

% nexttile
% s = pcolor(X,Y,alpha);
% s.EdgeColor = 'None';
% colorbar
% title('\alpha')

nexttile
s = pcolor(X,Y,beta);
s.EdgeColor = 'None';
colorbar
title('\beta')
drawnow

nexttile
s = pcolor(X,Y,real(rhs));
s.EdgeColor = 'None';
colorbar
title('rhs')
drawnow
 

% Constructing integral operator
gf = green(X - min(xs), Y - min(xs), b0 / a0, g0 / a0, false);
Gs = gf{1};
% fill in diagonal correction 
Gs_aug = [Gs, flip(Gs(1:end,2:end),2); ...
    flip(Gs(2:end,1:end)), flip(flip(Gs(2:end,2:end)),2)];
Gs_aug_hat = fft2(Gs_aug)*h*h;


% Solve with GMRES
start = tic;
mu = gmres(@(mu) lhs(mu,a0,bbar,Gs_aug_hat),rhs_vec,[],1e-12,200);
mu = reshape(mu, size(X));
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

N = n;
mu = reshape(mu, [N N]);
mu_aug = [mu, zeros(N,N-1); zeros(N-1,N), zeros(N-1)];
mu_aug_hat = fft2(mu_aug);
Gs_mu_aug = ifft2(Gs_aug_hat.*mu_aug_hat);
Gs_mu = Gs_mu_aug(1:N, 1:N);

% err = mu - bbar.*Gs_mu - k*bbar.*phiinc
 

% Evaluation and plotting
gf = green(X - min(xs), Y - min(xs), b0 / a0, g0 / a0, true);
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

phi_tot = phiinc + phi;
phi_n_tot = phi_n + k*phiinc;


figure(2);
tiledlayout(1,5)

nexttile
pc = pcolor(X,Y,real(mu));
pc.EdgeColor = 'none';
%clim([-1.5 1.5])
title('Re(\mu)')
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

ind = intersect(find(X == 4), find(Y == 4));
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

% Error of scattered part
% phi_n_sub = phi_n(ii-4:ii+4,jj-4:jj+4);
% first = a0*sum(bilap.*phi_n_sub,'all') ;
% second = -beta(ii,jj)*phi_n(ii,jj);
% third = g0*phi(ii,jj);
% rhs = k*bbar(ii,jj)*phiinc(ii,jj);
% err = abs(first + second + third - rhs) / max(abs(phi_n_tot(:)))

% Residual error of total solution 
phi_n_tot_sub = phi_n_tot(ii-4:ii+4,jj-4:jj+4);
first = a0*sum(bilap.*phi_n_tot_sub,'all') ;
second = -beta(ii,jj).*phi_n_tot(ii,jj);
third = g0*phi_tot(ii,jj);
err = abs(first + second + third) / max(abs(phi_n_tot(:)))

% Error of just the green's functions (?) Not really an error but a check
% Gs_sub = Gs(ii-4:ii+4,jj-4:jj+4);
% first = a0*sum(bilap.*Gs_sub,'all') ;
% second = -b0*Gs(ii,jj);
% third = g0*Gc(ii,jj);
% err = abs(first + second + third) 

% phiinc_sub = phiinc(ii-4:ii+4,jj-4:jj+4);
% first = k*sum(bilap.*phiinc_sub,'all') / h^4;
% second = -beta0*k*phiinc(ii,jj);
% third = gamma*phiinc(ii,jj);
% err = first + second + third

function v = lhs(mu,a0,bbar,Gs_aug_hat)
    N = sqrt(size(mu));
    N = N(1);
    mu = reshape(mu, [N N]);
    mu_aug = [mu, zeros(N,N-1); zeros(N-1,N), zeros(N-1)];
    mu_aug_hat = fft2(mu_aug);
    Gs_mu_aug = ifft2(Gs_aug_hat.*mu_aug_hat);
    Gs_mu = Gs_mu_aug(1:N, 1:N);
    v = a0*mu - bbar.*Gs_mu;
    v = v(:);
end