%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for plane wave 
% scattering of flexural-gravity waves with rolls
%
%%%%%

% clear 
% close all
% addpath(genpath('..'))

h = 20;

x1 = -4E3;
x2 = 10E3;

y1 = -7E3;
y2 = 7E3;

xs = x1:h:x2;
ys = y1:h:y2;
xl = 2*x1:h:2*x2;
yl = 2*y1:h:2*y2;
[~,n] = size(xs);
[X,Y] = meshgrid(xs,ys);
[XL,YL] = meshgrid(xl,yl);

freqs = [0.61;0.725];
ks = zeros(2,1);
iters = zeros(2,1);
phi_n_tots = zeros([size(X),2]);
mus = zeros([size(X),2]);

for ii = 1:numel(freqs)

    w = freqs(ii)

    % if w > 0.5
    %     h = 12.5;
    % end
    
    [coefs, H] = rolls(X,Y,0,40E2,-25E2,25E2,0.75,333.3,w); % remove gbar from coefs vector

    a0 = coefs{1}; 
    b0 = coefs{3}; 
    g0 = coefs{5}; 
    
    % Finding positive real roots
    [rts,ejs] = find_roots(b0 / a0, g0 / a0);
    k = rts((imag(rts) == 0) & (real(rts) > 0))
    ks(ii) = k;
    ejs = ejs/a0;
    
    % RHS (Incident field)
    k1 = k*cos(0);
    k2 = k*sin(0);
    phiinc = exp(1i*k1*X+1i*k2*Y);
    phininc = k*exp(1i*k1*X+1i*k2*Y);
    [rhs_vec, rhsp] = get_rhs_vec(coefs,k1,k2,phiinc);
    
    % Constructing integral operators
    src = [xl(ceil(end/2)); yl(ceil(end/2))];
    targ = [XL(:).'; YL(:).'];
    
    [inds,corrs] = get_correct(h,a0);
    kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs), h, inds,corrs);
    
    ind = find((XL == src(1)) & (YL == src(2)));
    sz = size(XL);
    
    kerns = gen_fft_kerns(kerns,sz,ind);
    
    evalkerns = {kerns{1}, kerns{4}};
    
    % Solve with GMRES
    [mu,flag,relres,iter,resvec] = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,[],1e-6,2000);
    iter
    iters(ii) = iter(2);
    mu = reshape(mu, size(X));
    
    [phi, phi_n] = sol_eval_fft(mu,evalkerns);

    
    phi_tot = phi + phiinc;
    phi_n_tot = phi_n + phininc;
    phi_n_tots(:,:,ii) = phi_n_tot;
    mus(:,:,ii) = mu;

    err = get_fin_diff_err(X,Y,mu,phi_n_tot,phi_tot,h,coefs,1760,1000)

end

%%

f= figure(1); clf 
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
f.Position = [70 303 726 371];

phi1 = phi_n_tots(:,:,1);
phi2 = phi_n_tots(:,:,2);

nexttile
pc = pcolor(X/1000,Y/1000,abs(phi1));
ylim([-6 6])
pc.EdgeColor = 'none';
pc.FaceColor = 'interp';
title('$k = 0.0235$ m$^{-1}$','FontWeight','normal','Interpreter','latex')
clim([0 max(abs(phi2(:)))])
drawnow
set(gca, 'FontSize',12)
xlabel('$x$ (km)','Interpreter','latex')
ylabel('$y$ (km)','Interpreter','latex')
axis square

nexttile
pc = pcolor(X/1000,Y/1000,abs(phi2));
ylim([-6 6])
pc.EdgeColor = 'none';
pc.FaceColor = 'interp';
title('$k = 0.0265$ m$^{-1}$','FontWeight','normal','Interpreter','latex')
colorbar
clim([0 max(abs(phi2(:)))])
drawnow
set(gca, 'FontSize',12)
xlabel('$x$ (km)','Interpreter','latex')
fontname(gcf, 'CMU Serif')
axis square

%% 
saveas(gcf,'rollfig3.fig','fig')
exportgraphics(gcf,'rollfig3.pdf','ContentType','image','Resolution',400)


return 
%% Figure generation for Jeremy

figure(4);

t = tiledlayout(1,4,"TileSpacing","tight");

X1 = X / 1000 + 2;
Y1 = Y / 1000 + 5;

ax1 = nexttile
s = pcolor(X1,Y1,H);
s.EdgeColor = 'None';
colormap(ax1,gray)
clim([min(H(:)) max(H(:))])
colorbar
title('H')


nexttile
pc = pcolor(X1,Y1,abs(mu));
clim([0 max(abs(mu(:)))])
pc.EdgeColor = 'none';
colorbar
title('|\rho|')

nexttile
pc = pcolor(X1,Y1,real(phi_tot));
pc.EdgeColor = 'none';
title('Re(\phi)')
colorbar

nexttile
pc = pcolor(X1,Y1,abs(phi_tot));
pc.EdgeColor = 'none';
title('|\phi|')
colorbar

%% Figure generation for Jeremy

figure(4);

t = tiledlayout('flow','TileSpacing','tight'); 

X1 = X / 1000 + 5;
Y1 = Y / 1000 + 5;

ax1 = nexttile
s = pcolor(X1,Y1,H*0+3);
s.EdgeColor = 'None';
colormap(ax1,gray)
clim([min(H(:)) max(H(:))])
colorbar
title('H')


nexttile
pc = pcolor(X1,Y1,H*0);
clim([0 max(abs(mu(:)))/3])
pc.EdgeColor = 'none';
colorbar
title('|\rho|')

nexttile([2 2]);
pc = pcolor(X1,Y1,H*0);
clim([0 max(abs(phi_n_tot(:)))*0.75])
pc.EdgeColor = 'none';
colorbar
title('|\phi_n|')


%% 

figure(5);
s = pcolor(X,Y,H);
s.EdgeColor = 'None';
colormap(gray)
%title('H')
axis off


%% 

figure(6);
pc = pcolor(X,Y,abs(mu));
clim([0 max(abs(mu(:)))/3])
pc.EdgeColor = 'none';
%title('|\rho|')
axis off

%% 

figure(7);
pc = pcolor(X,Y,abs(phi_n_tot));
clim([0 max(abs(phi_n_tot(:)))*0.75])
pc.EdgeColor = 'none';
%title('|\phi_n|')
axis off

%%

f= figure(1); 
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
f.Position = [70 303 726 371];

mu1 = mus(:,:,1);
mu2 = mus(:,:,2);

nexttile
pc = pcolor(X/1000,Y/1000,abs(mu1));
ylim([-6 6])
pc.EdgeColor = 'none';
pc.FaceColor = 'interp';
title('k = 0.0235','FontWeight','normal')
% clim([0 max(abs(mu2(:)))])
drawnow
set(gca, 'FontSize',12)
xlabel('x (km)')
ylabel('y (km)')
axis square

nexttile
pc = pcolor(X/1000,Y/1000,abs(mu2));
ylim([-6 6])
pc.EdgeColor = 'none';
pc.FaceColor = 'interp';
title('k = 0.0265','FontWeight','normal')
colorbar
% clim([0 max(abs(mu2(:)))])
drawnow
set(gca, 'FontSize',12)
xlabel('x (km)')
fontname(gcf, 'CMU Serif')
axis square
