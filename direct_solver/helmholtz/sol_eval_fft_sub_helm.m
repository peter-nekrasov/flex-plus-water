function phi_n = sol_eval_fft_sub_helm(mu,evalkerns,corr,h,dinds,iinds,jinds,X)

    Gs_aug_hat = evalkerns{1};

    Gs_corr = corr{1};

    N = size(X,1);

    dinds_aug = sub2ind([2*N-1,2*N-1], iinds, jinds);

    mu_aug = zeros([2*N-1,2*N-1]);
    mu_aug(dinds_aug) = mu;
    mu_aug_hat = fft2(mu_aug);

    mu0 = zeros(N^2,1);
    mu0(dinds) = mu;

    Gs_mu_aug = ifft2(Gs_aug_hat.*mu_aug_hat);
    Gs_mu = Gs_mu_aug(1:N,1:N);
    phi_n = Gs_mu(:) + Gs_corr*mu0*h*h;

    phi_n = reshape(phi_n,size(X));

end
