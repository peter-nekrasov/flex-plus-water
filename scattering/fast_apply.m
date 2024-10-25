function v = fast_apply(mu,kern_struct,V)

    alpha = V{1};
    alphax = V{2};
    alphay = V{3};
    alphaxx = V{4};
    alphaxy = V{5};
    alphayy = V{6};
    betabar = V{7};
    gammabar = V{8};
    nu = V{9};
    alphalap = alphaxx+alphayy;

    Gs_aug_hat = kern_struct{1};
    Gsxx_aug_hat = kern_struct{2};
    Gsxy_aug_hat = kern_struct{3};
    Gsyy_aug_hat = kern_struct{4};
    Gslapx_aug_hat = kern_struct{5};
    Gslapy_aug_hat = kern_struct{6};
    Gc_aug_hat = kern_struct{7};

    Gslap_aug_hat = Gsxx_aug_hat + Gsyy_aug_hat;

    N = sqrt(size(mu));
    N = N(1);
    
    mu = reshape(mu, [N N]);
    mu_aug = [mu, zeros(N,N-1); zeros(N-1,N), zeros(N-1)];
    mu_aug_hat = fft2(mu_aug);

    Gs_mu_aug = ifft2(Gs_aug_hat.*mu_aug_hat);
    Gs_mu = Gs_mu_aug(1:N, 1:N);

    Gsxx_mu_aug = ifft2(Gsxx_aug_hat.*mu_aug_hat);
    Gsxx_mu = Gsxx_mu_aug(1:N, 1:N);

    Gsxy_mu_aug = ifft2(Gsxy_aug_hat.*mu_aug_hat);
    Gsxy_mu = Gsxy_mu_aug(1:N, 1:N);

    Gsyy_mu_aug = ifft2(Gsyy_aug_hat.*mu_aug_hat);
    Gsyy_mu = Gsyy_mu_aug(1:N, 1:N);

    Gslap_mu_aug = ifft2(Gslap_aug_hat.*mu_aug_hat);
    Gslap_mu = Gslap_mu_aug(1:N, 1:N);

    Gslapx_mu_aug = ifft2(Gslapx_aug_hat.*mu_aug_hat);
    Gslapx_mu = Gslapx_mu_aug(1:N, 1:N);

    Gslapy_mu_aug = ifft2(Gslapy_aug_hat.*mu_aug_hat);
    Gslapy_mu = Gslapy_mu_aug(1:N, 1:N);

    Gc_mu_aug = ifft2(Gc_aug_hat.*mu_aug_hat);
    Gc_mu = Gc_mu_aug(1:N, 1:N);

    % write function to apply integral operators above

    v = alpha.*mu + alphalap.*Gslap_mu + 2*alphax.*Gslapx_mu + ...
        + 2*alphay.*Gslapy_mu + (1-nu)*(2*alphaxy.*Gsxy_mu - ...
        alphaxx.*Gsyy_mu - alphayy.*Gsxx_mu) - betabar.*Gs_mu + gammabar.*Gc_mu;
    v = v(:);

end