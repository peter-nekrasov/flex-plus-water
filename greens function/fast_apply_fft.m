function v = fast_apply_fft(mu,kern_struct,V)

    a0 = V{1};
    abar = V{2};
    b0 = V{3};
    bbar = V{4};
    g0 = V{5};
    gbar = V{6};
    alphax = V{7};
    alphay = V{8};

    % alphaxx = V{9};
    % alphaxy = V{10};
    % alphayy = V{11};
    % nu = V{12};

    %alphalap = alphaxx+alphayy;

    Gs_aug_hat = kern_struct{1};
    Gs_hess = kern_struct{2};
    Gs_gradlap = kern_struct{3};
    Gphi_aug_hat = kern_struct{4};

    Gsxx_aug_hat = Gs_hess(:,:,1);
    Gsxy_aug_hat = Gs_hess(:,:,2);
    Gsyy_aug_hat = Gs_hess(:,:,3);

    Gslapx_aug_hat = Gs_gradlap(:,:,1);
    Gslapy_aug_hat = Gs_gradlap(:,:,2);

    Gslap_aug_hat = Gsxx_aug_hat + Gsyy_aug_hat;

    N = sqrt(size(mu));
    N = N(1);
    
    mu = reshape(mu, [N N]);
    mu_aug = [mu, zeros(N,N-1); zeros(N-1,N), zeros(N-1)];
    mu_aug_hat = fft2(mu_aug);

    Gs_mu_aug = ifft2(Gs_aug_hat.*mu_aug_hat);
    Gs_mu = Gs_mu_aug(1:N, 1:N);

    Gphi_mu_aug = ifft2(Gphi_aug_hat.*mu_aug_hat);
    Gphi_mu = Gphi_mu_aug(1:N, 1:N);

    % Gsxx_mu_aug = ifft2(Gsxx_aug_hat.*mu_aug_hat);
    % Gsxx_mu = Gsxx_mu_aug(1:N, 1:N);
    % 
    % Gsxy_mu_aug = ifft2(Gsxy_aug_hat.*mu_aug_hat);
    % Gsxy_mu = Gsxy_mu_aug(1:N, 1:N);
    % 
    % Gsyy_mu_aug = ifft2(Gsyy_aug_hat.*mu_aug_hat);
    % Gsyy_mu = Gsyy_mu_aug(1:N, 1:N);
    % 
    % Gslap_mu_aug = ifft2(Gslap_aug_hat.*mu_aug_hat);
    % Gslap_mu = Gslap_mu_aug(1:N, 1:N);

    Gslapx_mu_aug = ifft2(Gslapx_aug_hat.*mu_aug_hat);
    Gslapx_mu = Gslapx_mu_aug(1:N, 1:N);

    Gslapy_mu_aug = ifft2(Gslapy_aug_hat.*mu_aug_hat);
    Gslapy_mu = Gslapy_mu_aug(1:N, 1:N);

    % Gc_mu_aug = ifft2(Gc_aug_hat.*mu_aug_hat);
    % Gc_mu = Gc_mu_aug(1:N, 1:N);


    % v = alpha.*mu + alphalap.*Gslap_mu + 2*alphax.*Gslapx_mu + ...
    %     + 2*alphay.*Gslapy_mu + (1-nu)*(2*alphaxy.*Gsxy_mu - ...
    %     alphaxx.*Gsyy_mu - alphayy.*Gsxx_mu) - betabar.*Gs_mu + gammabar.*Gc_mu;
    % v = v(:);

    v = (a0 + abar)./a0.*mu + alphax.*Gslapx_mu + ...
        + alphay.*Gslapy_mu - 0.5*bbar.*Gs_mu + 0.5*gbar.*Gphi_mu ;
    v = v(:);

end