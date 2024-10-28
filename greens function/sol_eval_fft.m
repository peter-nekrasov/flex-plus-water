function [phi, phi_n] = sol_eval_fft(mu,evalkerns)

    [n,~] = size(mu);

    Gs_aug_hat = evalkerns{1};
    Gc_aug_hat = evalkerns{2};
        
    % Evaluation and plotting
    mu_aug = [mu, zeros(n,n-1); zeros(n-1,n), zeros(n-1)];
    mu_aug_hat = fft2(mu_aug);
    phi_aug = ifft2(Gc_aug_hat.*mu_aug_hat);
    phi_n_aug = ifft2(Gs_aug_hat.*mu_aug_hat);
    phi = phi_aug(1:n,1:n);
    phi_n = phi_n_aug(1:n,1:n);
    
end
