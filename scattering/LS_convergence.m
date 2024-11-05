%%%%%
%
% Testing convergence of the adjointed Lippman-Schwinger equation 
% using finite differences
%
%%%%%

clear 
close all
addpath(genpath('..'))

L = 50;
hs = [5 2.5 1 0.5 0.25 0.2 0.125 0.1 0.125/2 5/90];
errs = hs*0;

for ii = 1:numel(hs)

h = hs(ii);

xs = -L:h:L;
xl = -2*L:h:2*L;
[~,n] = size(xs);
[X,Y] = meshgrid(xs);
[XL,YL] = meshgrid(xl);

coefs = bump(X,Y,4,6); % remove gbar from coefs vector

a0 = coefs{1}; 
b0 = coefs{3}; 
g0 = coefs{5}; 

% Finding positive real roots
[rts,ejs] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0));
ejs = ejs/a0;

src = [0;0];
targ = [XL(:).'; YL(:).'];

% RHS (Incident field)
k1 = 2*sqrt(2)*k/3;
k2 = k/3;
phiinc = exp(1i*k1*X+1i*k2*Y);
[rhs_vec, rhs] = get_rhs_vec(coefs,k1,k2,phiinc);

% figure(1);
% tiledlayout(1,3);
%
% nexttile
% s = pcolor(X,Y,coefs{1} + coefs{2});
% s.EdgeColor = 'None';
% colorbar
% title('\alpha')
% drawnow
% 
% nexttile
% s = pcolor(X,Y,coefs{2} + coefs{3});
% s.EdgeColor = 'None';
% colorbar
% title('\beta')
% drawnow
% 
% nexttile
% s = pcolor(X,Y,real(rhs));
% s.EdgeColor = 'None';
% colorbar
% title('rhs')
% drawnow


% Constructing integral operators
[inds,corrs] = get_correct(rts,ejs,h,a0);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs), inds,corrs);

ind = find((XL == 0) & (YL ==0));
sz = size(XL);

kerns = gen_fft_kerns(kerns,sz,ind);

evalkerns = {kerns{1}, kerns{4}};

% Solve with GMRES
start = tic;
mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,[],1e-12,200);
mu = reshape(mu, size(X));
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

[phi, phi_n] = sol_eval_fft(mu,evalkerns);

phi_tot = phi + phiinc;
phi_n_tot = phi_n + k*phiinc;


% figure(2);
% tiledlayout(1,5)
% 
% nexttile
% pc = pcolor(X,Y,real(mu));
% pc.EdgeColor = 'none';
% title('Re(\mu)')
% colorbar
% 
% nexttile
% pc = pcolor(X,Y,real(phi_tot));
% pc.EdgeColor = 'none';
% title('Re(\phi)')
% colorbar
% 
% nexttile
% pc = pcolor(X,Y,abs(phi_tot));
% pc.EdgeColor = 'none';
% title('|\phi|')
% colorbar
% 
% nexttile
% pc = pcolor(X,Y,real(phi_n_tot));
% pc.EdgeColor = 'none';
% title('real(\phi_n)')
% colorbar
% 
% nexttile
% pc = pcolor(X,Y,abs(phi_n_tot));
% pc.EdgeColor = 'none';
% title('|\phi_n|')
% colorbar
       
% Calculate error with finite difference
errs(ii) = get_fin_diff_err(X,Y,mu,phi_n_tot,phi_tot,h,coefs)

end

%% Plotting
 
figure(1);
t = tiledlayout(1,2);
title(t, ['\alpha =  ', num2str(a0),', \beta = ',num2str(b0), ', \gamma = ',num2str(g0)]);

nexttile
s = pcolor(X,Y,coefs{1} + coefs{2});
s.EdgeColor = 'None';
colorbar
title('\alpha')
hold on
scatter(5,5,75,"MarkerEdgeColor","k", ...
    "MarkerFaceColor","k")


nexttile
N = 100./hs(1:end-1);
loglog(N,errs(1:end-1),'o-')
hold on
loglog(N,10000000*N.^(-6),'--')
xlim([0.5*min(N), 2*max(N)])
ylim([0.01*min(errs(1:end-1)) 100*max(errs(1:end-1))])
xlabel('N')
ylabel('Relative error')
title('Residual')
legend('errors','N^{-6}')
