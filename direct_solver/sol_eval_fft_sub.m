function [phi, phi_n] = sol_eval_fft_sub(mu,evalkerns,corr,h,dinds,iinds,jinds,X)

    Gs_aug_hat = evalkerns{1};
    Gphi_aug_hat = evalkerns{2};

    Gs_corr = corr{1};
    Gphi_corr = corr{2};

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

    Gphi_mu_aug = ifft2(Gphi_aug_hat.*mu_aug_hat);
    Gphi_mu = Gphi_mu_aug(1:N,1:N);
    phi = Gphi_mu(:) + Gphi_corr*mu0*h*h;
        
    phi = phi/2;
    phi_n = phi_n/2;

    phi = reshape(phi,size(X));
    phi_n = reshape(phi_n,size(X));

end
