function v = fast_apply_fft_sub_helm(mu,kern_struct,zk,coefs,corr,h,dinds,iinds,jinds,X)

    V = coefs{1}(dinds);

    Gs_aug_hat = kern_struct{1};
    Gs_corr = corr{1}(dinds,dinds);

    N = size(X,1);

    dinds_aug = sub2ind([2*N-1,2*N-1], iinds, jinds);

    mu_aug = zeros([2*N-1,2*N-1]);
    mu_aug(dinds_aug) = mu;
    mu_aug_hat = fft2(mu_aug);

    Gs_mu_aug = ifft2(Gs_aug_hat.*mu_aug_hat);
    Gs_mu = Gs_mu_aug(dinds_aug);
    Gs_mu = Gs_mu(:) + Gs_corr*mu*h*h;

    v = mu + zk^2*V.*Gs_mu;

end