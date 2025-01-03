clear 
close all
addpath(genpath('..'))

h = 0.02;
L = 7;

xs = -L/2:h:L/2;
xl = -L:h:L;

ys = -L/2:h:L/2;
yl = -L:h:L;

% [~,~,coefs,Hs] = spiral(-L/2,L/2,-L/2,L/2,h,-0.7,50);
%[~,~,coefs,Hs] = spiral(-L/2,L/2,-L/2,L/2,h,-0.75,100);
%[~,~,coefs,Hs] = spiral(-L/2,L/2,-L/2,L/2,h,-0.8,15);
%[~,~,coefs,Hs] = logspiral(-L/2,L/2,-L/2,L/2,h,-0.7,300);
[~,~,coefs,Hs] = logspiral(-L/2,L/2,-L/2,L/2,h,-0.6,750);


[X,Y] = meshgrid(xs,ys);
[XL,YL] = meshgrid(xl,yl);

H = Hs{1};

a0 = coefs{1}; 
b0 = coefs{3}; 
g0 = coefs{5}; 

% Finding positive real roots
[rts,ejs] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0))
ejs = ejs/a0;



%%
src = [0;0];
targ = [XL(:).'; YL(:).'];

% RHS (Incident field)
k1 = k*cos(pi/3);
k2 = k*sin(pi/3);
phiinc = exp(1i*k1*X+1i*k2*Y);
[rhs_vec, rhs] = get_rhs_vec(coefs,k1,k2,phiinc);

nexttile
s = pcolor(X,Y,H);
s.EdgeColor = 'None';
colorbar
title('H')
drawnow

nexttile
s = pcolor(X,Y,(coefs{1} + coefs{2}));
s.EdgeColor = 'None';
colorbar
title('\alpha')
drawnow

nexttile
s = pcolor(X,Y,(coefs{2} + coefs{3}));
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
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h, inds,corrs);

ind = find((XL == 0) & (YL ==0));
sz = size(XL);

kerns = gen_fft_kerns(kerns,sz,ind);

evalkerns = {kerns{1}, kerns{4}};

% mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,30,1e-7,1000);
mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,30,1e-7,2000);
mu = reshape(mu, size(X));

[phi, phi_n] = sol_eval_fft(mu,evalkerns);


phiinc = exp(1i*k1*X+1i*k2*Y);
phininc = k*exp(1i*k1*X+1i*k2*Y);

phi_tot = phi + phiinc;
phi_n_tot = phi_n + phininc;

load gong.mat
sound(y)

%% 

figure(2);
tiledlayout(1,3,'TileSpacing','tight')
ax1 = nexttile;
pl1 = pcolor(X+4,Y+4,H);
pl1.EdgeColor = 'none';
title('Thickness (m)','FontWeight','normal')
c = [0.3:0.01:0.9 ; 0.3:0.01:0.9; 0.3:0.01:0.9 ].';
% xlim([0 9])
% ylim([0 9])
ylabel('y (m)')
xlabel('x (m)')
colormap(ax1,c)
colorbar
axis square
set(gca, 'FontSize',12)


nexttile
pl2 = pcolor(X+4,Y+4,real(phi_n_tot/k));
pl2.EdgeColor = 'none';
title('Re(\phi_n)','FontWeight','normal')
colorbar
% xlim([0 9])
% ylim([0 9])
xlabel('x (m)')
%clim([-3 3])
axis square
set(gca, 'FontSize',12)



nexttile
pl2 = pcolor(X+4,Y+4,abs(phi_n_tot/k));
pl2.EdgeColor = 'none';
title('|\phi_n|','FontWeight','normal')
colorbar
% xlim([0 9])
% ylim([0 9])
xlabel('x (m)')
%clim([-3 3])
axis square

set(gca, 'FontSize',12)
xlabel('x (m)')

fontname(gcf, 'CMU Serif')

return

%%

saveas(gcf,'spiral2.fig','fig')
exportgraphics(gcf,'spiral2.pdf','ContentType','vector','Resolution',1200)