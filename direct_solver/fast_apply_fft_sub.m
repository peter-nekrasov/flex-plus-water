function v = fast_apply_fft_sub(mu,kern_struct,V,corr,h,dinds,iinds,jinds,X)

    a0 = V{1};
    abar = V{2}(dinds);
    b0 = V{3};
    bbar = V{4}(dinds);
    g0 = V{5};
    alphax = V{7}(dinds);
    alphay = V{8}(dinds);
    alphaxx = V{9}(dinds);
    alphaxy = V{10}(dinds);
    alphayy = V{11}(dinds);
    nu = V{end};
    alphalap = alphaxx + alphayy;

    Gs_aug_hat = kern_struct{1};
    Gs_hess = kern_struct{2};
    Gs_gradlap = kern_struct{3};
    Gphi_aug_hat = kern_struct{4};

    Gsxx_aug_hat = Gs_hess(:,:,1);
    Gsxy_aug_hat = Gs_hess(:,:,2);
    Gsyy_aug_hat = Gs_hess(:,:,3);

    Gslapx_aug_hat = Gs_gradlap(:,:,1);
    Gslapy_aug_hat = Gs_gradlap(:,:,2);

    Gs_corr = corr{1}(dinds,dinds);
    Gs_hess_corr = corr{2};
    Gs_gradlap_corr = corr{3};
    Gphi_corr = corr{4}(dinds,dinds);

    Gs_xx_corr = Gs_hess_corr{1}(dinds,dinds);
    Gs_xy_corr = Gs_hess_corr{2}(dinds,dinds);
    Gs_yy_corr = Gs_hess_corr{3}(dinds,dinds);

    Gslapx_corr = Gs_gradlap_corr{1}(dinds,dinds);
    Gslapy_corr = Gs_gradlap_corr{2}(dinds,dinds);

    Gslap_corr = Gs_xx_corr + Gs_yy_corr;
    Gslap_aug_hat = Gsxx_aug_hat + Gsyy_aug_hat;

    N = size(X,1);

    dinds_aug = sub2ind([2*N-1,2*N-1], iinds, jinds);

    mu_aug = zeros([2*N-1,2*N-1]);
    mu_aug(dinds_aug) = mu;
    mu_aug_hat = fft2(mu_aug);

    Gs_mu_aug = ifft2(Gs_aug_hat.*mu_aug_hat);
    Gs_mu = Gs_mu_aug(dinds_aug);
    Gs_mu = Gs_mu(:) + Gs_corr*mu*h*h;

    Gphi_mu_aug = ifft2(Gphi_aug_hat.*mu_aug_hat);
    Gphi_mu = Gphi_mu_aug(dinds_aug);
    Gphi_mu = Gphi_mu(:) + Gphi_corr*mu*h*h;

    Gsxx_mu_aug = ifft2(Gsxx_aug_hat.*mu_aug_hat);
    Gsxx_mu = Gsxx_mu_aug(dinds_aug);
    Gsxx_mu = Gsxx_mu(:) + Gs_xx_corr*mu*h*h;

    Gsxy_mu_aug = ifft2(Gsxy_aug_hat.*mu_aug_hat);
    Gsxy_mu = Gsxy_mu_aug(dinds_aug);
    Gsxy_mu = Gsxy_mu(:) + Gs_xy_corr*mu*h*h;

    Gsyy_mu_aug = ifft2(Gsyy_aug_hat.*mu_aug_hat);
    Gsyy_mu = Gsyy_mu_aug(dinds_aug);
    Gsyy_mu = Gsyy_mu(:) + Gs_yy_corr*mu*h*h;

    Gslap_mu_aug = ifft2(Gslap_aug_hat.*mu_aug_hat);
    Gslap_mu = Gslap_mu_aug(dinds_aug);
    Gslap_mu = Gslap_mu(:) + Gslap_corr*mu*h*h;

    Gslapx_mu_aug = ifft2(Gslapx_aug_hat.*mu_aug_hat);
    Gslapx_mu = Gslapx_mu_aug(dinds_aug); 
    Gslapx_mu = Gslapx_mu(:) + Gslapx_corr*mu*h*h;

    Gslapy_mu_aug = ifft2(Gslapy_aug_hat.*mu_aug_hat);
    Gslapy_mu = Gslapy_mu_aug(dinds_aug);
    Gslapy_mu = Gslapy_mu(:) + Gslapy_corr*mu*h*h;

    v = (a0 + abar)./a0.*mu + alphax.*Gslapx_mu + alphay.*Gslapy_mu ...
        + 0.5.*alphalap.*Gslap_mu + ...
        + 0.5*(1-nu)*(2*alphaxy.*Gsxy_mu-alphayy.*Gsxx_mu-alphaxx.*Gsyy_mu) ...
        - 0.5*(a0*bbar-abar*b0)./a0.*Gs_mu - 0.5*abar.*g0./a0.*Gphi_mu ;
    v = v(:);

end