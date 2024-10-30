% Parameters
h = 0.2;
xs = -50:h:50; %-5000:h:5000;
[~,N] = size(xs);
[X,Y] = meshgrid(xs);
% H0 = 20;
a0 = 1; %6.410256410256411e+08*H0^3;
b0 = 3; %(917*H0*w^2 - 9800) ;
g0 = -1; %- 1000*w^2 ;

dens = @(x,y) exp(-((x).^2 + (y).^2)/(2*(4*k)^2));
testf = @(x,y) greenval(x,y,b0,g0,true).*dens(x,y);
val = integral2(testf,-50,50,-50,50);

% Finding positive real roots
[rts,~] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0));

% Constructing integral operator
gf = green(X - min(xs), Y - min(xs), b0 / a0, g0 / a0, true);
Gs = gf{1};
% fill in diagonal correction 
Gs_aug = [Gs, flip(Gs(1:end,2:end),2); ...
    flip(Gs(2:end,1:end)), flip(flip(Gs(2:end,2:end)),2)];
Gs_aug_hat = fft2(Gs_aug)*h*h;

% Creating density
mu = dens(X,Y);

% Integrating
mu_aug = [mu, zeros(N,N-1); zeros(N-1,N), zeros(N-1)];
mu_aug_hat = fft2(mu_aug);
Gs_mu_aug = ifft2(Gs_aug_hat.*mu_aug_hat);
Gs_mu = Gs_mu_aug(1:N, 1:N);

ind = find(X == 0 & Y == 0);
err = abs(Gs_mu(ind) - val) / abs(val)


