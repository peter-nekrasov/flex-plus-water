%%%%%
%
% Continuous scattering with variable beta 
%
%%%%%

addpath(genpath('..'))

L = 50;
hs = 5./2.^(1:6);
errs = hs*0;

% H0 = 20;
% w = 2;
a0 = 1; %6.410256410256411e+08*H0^3;
b0 = 3; %(917*H0*w^2 - 9800) ;
g0 = -1; %- 1000*w^2 ;

% Finding positive real roots
[rts,~] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0));

for ii = 1:numel(hs)

    % Parameters
    h = hs(ii);

    xs = -L:h:L;
    xl = -2*L:h:2*L;
    [~,n] = size(xs);
    [X,Y] = meshgrid(xs);
    [XL,YL] = meshgrid(xl);

    bbar = -2*X.*exp(-(X.^2 + Y.^2)/(2*(4*k)^2));
    beta = b0 + bbar;
    
    coefs = {a0,bbar};
    
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
     
    % Finding diagonal entry and shifting kernels
    %[zi,zj] = ind2sub(size(X),ind);
    
    % Constructing integral operators
    Gs = green(XL,YL, b0 / a0, g0 / a0, false);
    Gc = green(XL,YL, b0 / a0, g0 / a0, true);
    %Gs = green(X - min(xs),Y - min(xs), b0 / a0, g0 / a0, false);
    %Gc = green(X - min(xs),Y - min(xs), b0 / a0, g0 / a0, true);

    kerns = [Gs, Gc{1}];
    ind = find((XL == 0) & (YL ==0));
    kerns = proc_kern(kerns,h,ind,b0/a0,g0/a0);
    
    Gs_aug_hat = kerns{1};
    Gc_aug_hat = kerns{7};
    
    % Solve with GMRES
    start = tic;
    mu = gmres(@(mu) lhs(mu,kerns,coefs),rhs_vec,[],1e-12,200);
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
           
    % Calculate error with finite difference
    errs(ii) = get_fin_diff_err(X,Y,mu,phi_n_tot,phi_tot,a0,beta,g0,h);

end

%% Plotting errors 

figure(3)
loglog(hs,errs,'x-');

hold on
loglog(hs, 0.01*hs.^8, '--')

