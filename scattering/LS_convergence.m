%%%%%
%
% Testing convergence of the adjointed Lippman-Schwinger equation 
% using finite differences
%
%%%%%

clear 
close all
addpath(genpath('..'))

L = 4000;
hs = [200 100 50 25 20 10 8 5 4 10/3];
errsinf = hs*0;
errs2 = hs*0;

for ii = 1:numel(hs)

h = hs(ii);

xs = -L:h:L;
xl = -2*L:h:2*L;
[~,n] = size(xs);
[X,Y] = meshgrid(xs);
[XL,YL] = meshgrid(xl);

[coefs, H] = bump2(X,Y,5,200,1); % remove gbar from coefs vector

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

figure(1);
tiledlayout(1,3);

nexttile
s = pcolor(X,Y,coefs{1} + coefs{2});
s.EdgeColor = 'None';
colorbar
title('\alpha')
drawnow

nexttile
s = pcolor(X,Y,coefs{2} + coefs{3});
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


% Constructing integral operators
[inds,corrs] = get_correct(h,a0);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs), h, inds,corrs);

ind = find((XL == 0) & (YL ==0));
sz = size(XL);

kerns = gen_fft_kerns(kerns,sz,ind);

evalkerns = {kerns{1}, kerns{4}};

% Solve with GMRES
start = tic;
mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,30,1e-14,200);
mu = reshape(mu, size(X));
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

[phi, phi_n] = sol_eval_fft(mu,evalkerns);

phi_tot = phi + phiinc;
phi_n_tot = phi_n + k*phiinc;


figure(2);
tiledlayout(1,5)

nexttile
pc = pcolor(X,Y,real(mu));
pc.EdgeColor = 'none';
title('Re(\mu)')
colorbar

nexttile
pc = pcolor(X,Y,real(phi_tot));
pc.EdgeColor = 'none';
title('Re(\phi)')
colorbar

nexttile
pc = pcolor(X,Y,abs(phi_tot));
pc.EdgeColor = 'none';
title('|\phi|')
colorbar

nexttile
pc = pcolor(X,Y,real(phi_n_tot));
pc.EdgeColor = 'none';
title('real(\phi_n)')
colorbar

nexttile
pc = pcolor(X,Y,abs(phi_n_tot));
pc.EdgeColor = 'none';
title('|\phi_n|')
colorbar
       
[Xerr,Yerr] = meshgrid(-L/2:400:L/2);

% Calculate error with finite difference
err = get_fin_diff_err3(X,Y,mu,phi_n_tot,phi_tot,h,coefs,Xerr(:),Yerr(:));

errsinf(ii) = max(err);
errs2(ii) = mean(err);

end

%% Plotting
 
figure(4);


nexttile
loglog(hs,errs2,'o-')
hold on
loglog(hs,errsinf,'o-')
hold on

loglog(hs(3:end-2),4E-17*hs(3:end-2).^6,'k--')
hold on


xlim([0.5*min(hs), 2*max(hs)])
ylim([0.01*min(errs2(1:end-1)) 100*max(errs2(1:end-1))])
xlabel('N')
ylabel('Relative error')
title('Residual')
legend('$\ell_2$ - fin diff','$\ell_\infty$ - fin diff','$h^6$','Interpreter','latex')
