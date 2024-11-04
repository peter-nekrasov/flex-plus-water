%%%%%
%
% Solving the Adjointed-Lippman Schwinger equation
%
%%%%%

addpath('../greens function/')
addpath('geometry/')

% Parameters
h = 0.5;
xs = -50:h:50;
[~,n] = size(xs);
[X,Y] = meshgrid(xs);

% Physical constants
alpha0 = 10;
beta0 = 5;
gamma0 = -2;
nu = 0.3;

% Calculating perturbed coefficients
alpha = gaussian(X,Y);
gammabar = - alpha{1}*gamma0; % check these are right
betabar = - gamma0*alpha0*real((alpha{1})^(1/3)) + alpha{1}*beta0; % check these are right
alpha{1} = alpha{1} + alpha0;
gamma = gamma0 + gammabar;
beta = beta0 + betabar;
coefs = [alpha {betabar,gammabar,nu}];

figure(1);
tiledlayout(1,3);

nexttile
s = pcolor(X,Y,alpha{1});
s.EdgeColor = 'None';
colorbar
title('\alpha')

nexttile
s = pcolor(X,Y,beta);
s.EdgeColor = 'None';
colorbar
title('\beta')

nexttile
s = pcolor(X,Y,gamma);
s.EdgeColor = 'None';
colorbar
title('\gamma')
drawnow

figure(2);
tiledlayout(2,3);

nexttile
s = pcolor(X,Y,alpha{2});
s.EdgeColor = 'None';
colorbar
title('\partial_x \alpha')

nexttile
s = pcolor(X,Y,alpha{3});
s.EdgeColor = 'None';
colorbar
title('\partial_y \alpha')

nexttile
s = pcolor(X,Y,alpha{4});
s.EdgeColor = 'None';
colorbar
title('\partial_{xx} \alpha')

nexttile
s = pcolor(X,Y,alpha{5});
s.EdgeColor = 'None';
colorbar
title('\partial_{xy} \alpha')

nexttile
s = pcolor(X,Y,alpha{6});
s.EdgeColor = 'None';
colorbar
title('\partial_{yy} \alpha')

nexttile
s = pcolor(X,Y,alpha{4}+alpha{6});
s.EdgeColor = 'None';
colorbar
title('\Delta \alpha')
drawnow
 
% Finding positive real roots
[rts,~] = find_roots(beta0,gamma0);
k = rts((imag(rts) == 0) & (real(rts) > 0));

% Finding diagonal entry and shifting kernels
ind = find((X == 0) & (Y==0));
[zi,zj] = ind2sub(size(X),ind);

Xshift = circshift(X,[zi,zj]);
Yshift = circshift(Y,[zi,zj]);

% Constructing integral operators
Gs = green(Xshift,Yshift, beta0, gamma0, false);
Gc = green(Xshift,Yshift, beta0, gamma0, true);
kerns = [Gs, Gc{1}];
kerns = proc_kern(kerns,h); % the singularity has to be at (1,1)

% RHS (Incident field)
phiinc = exp(1i*k*X);
alphalap = alpha{4} + alpha{6};
alphax = alpha{2};
alphayy = alpha{6};
rhs = (-k*(k^2*alphalap+2*k^3*alphax+(1-nu)*(-k^2*alphayy)) + k*betabar - gammabar).*phiinc;
rhs_vec = rhs(:);

figure(3)
s = pcolor(X,Y,real(rhs));
s.EdgeColor = 'None';
colorbar
title('RHS vector')
drawnow

% Solve with GMRES
start = tic;
mu = gmres(@(mu) lhs(mu,kerns,coefs),rhs_vec,[],1e-12,200);
mu = reshape(mu, size(X));
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

% Evaluation and plotting
mu_aug = [mu, zeros(n,n-1); zeros(n-1,n), zeros(n-1)];
mu_aug_hat = fft2(mu_aug);
phi_aug = ifft2(kerns{7}.*mu_aug_hat);
phi_n_aug = ifft2(kerns{1}.*mu_aug_hat);
phi = phi_aug(1:n,1:n);
phi_n = phi_n_aug(1:n,1:n);

phi_tot = phiinc+phi;
phi_n_tot = phi_n + k*phiinc;

figure(4)
tiledlayout(1,4)

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

% ind = intersect(find(X == 16), find(Y == 0));
% [ii, jj] = ind2sub(size(X),ind);
% 
% % finite difference stencils
% d1 = zeros(9,1);
% d1(1) = 7/240;	
% d1(2) = -2/5;
% d1(3) = 169/60;
% d1(4) = -122/15;
% d1(5) = 91/8;
% d1(6) = -122/15;
% d1(7) = 169/60;
% d1(8) = -2/5;
% d1(9) = 7/240;
% 
% d2 = zeros(9, 1);
% d2(1) = -1/560;
% d2(2) = 8/315;
% d2(3) = -1/5;
% d2(4) = 8/5;
% d2(5) = -205/72;
% d2(6) = 8/5;
% d2(7) = -1/5;
% d2(8) = 8/315;
% d2(9) = -1/560;
% 
% bilap = zeros(9,9);
% bilap(:,5) = d1;
% bilap(5,:) = bilap(5,:) + d1.';
% bilap = bilap + 2*(d2*d2.');
% bilap = bilap / h^4;
% 
% phi_n_tot_sub = phi_n_tot(ii-4:ii+4,jj-4:jj+4);
% first = sum(bilap.*phi_n_tot_sub,'all') ;
% second = -beta(ii,jj)*phi_n_tot(ii,jj);
% third = gamma(ii,jj)*phi_tot(ii,jj);
% err = abs(first + second + third) / max(abs(phi_n(:)))

%phiinc_sub = phiinc(ii-4:ii+4,jj-4:jj+4);
%first = k*sum(bilap.*phiinc_sub,'all') / h^4;
%second = -beta0*k*phiinc(ii,jj);
%third = gamma*phiinc(ii,jj);
%err = first + second + third

function v = lhs(mu,kern_struct,V)

    v = fast_apply(mu,kern_struct,V);

end