clear 
close all
addpath(genpath('..'))

%h = 2;
%L = 1100;

h = 5;
L = 1700;

xs = -L/2:h:L/2;
xl = -L:h:L;

ys = -L/2:h:L/2;
yl = -L:h:L;

[~,~,coefs,Hs] = spiral(-L/2,L/2,-L/2,L/2,h,0.66,2.9);
% [~,~,coefs,Hs] = spiral(-L/2,L/2,-L/2,L/2,h,-0.65,3);
%[~,~,coefs,Hs] = spiral(-L/2,L/2,-L/2,L/2,h,-0.67,4.8);
% [~,~,coefs,Hs] = spiral(-L/2,L/2,-L/2,L/2,h,-0.6,5.5);
%[~,~,coefs,Hs] = spiral(-L/2,L/2,-L/2,L/2,h,-0.75,100);
%[~,~,coefs,Hs] = spiral(-L/2,L/2,-L/2,L/2,h,-0.8,15);

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
k1 = k*cos(2*pi/3);
k2 = k*sin(2*pi/3);
k1 = k*cos(2*pi/3);
k2 = k*sin(2*pi/3);
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
mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,30,1e-5,500);
mu = reshape(mu, size(X));

[phi, phi_n] = sol_eval_fft(mu,evalkerns);


phiinc = exp(1i*k1*X+1i*k2*Y);
phininc = k*exp(1i*k1*X+1i*k2*Y);

phi_tot = phi + phiinc;
phi_n_tot = phi_n + phininc;

load gong.mat
sound(y)


%%

f = figure(2);
f.Units = 'points';
f.InnerPosition = [584 281 700 634];

figure(2);
tiledlayout(2,2,'TileSpacing','tight')
ax1 = nexttile;
pl1 = pcolor(X/1000+0.67,Y/1000+0.67,H,'FaceColor','interp');
pl1.EdgeColor = 'none';
title('Thickness (m)','FontWeight','normal')
c1 = [0.7*ones(1,71); 0:0.01:0.7; 0:0.01:0.7  ].' / 0.7;
c2 = [ 2:-0.01:0; 2:-0.01:0; 2*ones(1,201);].'/2;
c = [c1; c2];
xlim([0 1.5])
ylim([0 1.5])
ylabel('y (km)')
clim([0.3 3])
colormap(ax1,c)
colorbar
axis square
set(gca, 'FontSize',12)


nexttile
pl2 = pcolor(X/1000+0.67,Y/1000+0.67,abs(mu),'FaceColor','interp');
pl2.EdgeColor = 'none';
title('|\mu|','FontWeight','normal')
colorbar
xlim([0 1.5])
ylim([0 1.5])
%clim([-3 3])
axis square
set(gca, 'FontSize',12)

annotation('arrow',[0.905 0.855],[0.111 0.2])
annotation('arrow',[0.905-0.442 0.855-0.442],[0.111 0.2])


nexttile
pl2 = pcolor(X/1000+0.67,Y/1000+0.67,real(phi_n_tot),'FaceColor','interp');
pl2.EdgeColor = 'none';
title('\Re(\phi_z)','FontWeight','normal')
colorbar
xlim([0 1.5])
ylim([0 1.5])
ylabel('y (km)')
xlabel('x (km)')
clim([0.7*min(real(phi_n_tot(:))) 0.7*max(real(phi_n_tot(:)))])
axis square
set(gca, 'FontSize',12)
xlabel('x (km)')


nexttile
pl2 = pcolor(X/1000+0.67,Y/1000+0.67,abs(phi_n_tot),'FaceColor','interp');
pl2.EdgeColor = 'none';
title('|\phi_z|','FontWeight','normal')
colorbar
xlim([0 1.5])
ylim([0 1.5])
xlabel('x (km)')
%clim([-3 3])
clim([0 0.9*max(abs(phi_n_tot(:)))])
axis square
set(gca, 'FontSize',12)
xlabel('x (km)')

fontname(gcf, 'CMU Serif')

return



%% 
% f = figure(2);
% f.Units = 'points';
% f.InnerPosition = [584 281 700 634];
% 
% figure(2);
% tiledlayout(2,2,'TileSpacing','tight')
% ax1 = nexttile;
% pl1 = pcolor(X/1000+0.45,Y/1000+0.45,H,'FaceColor','interp');
% pl1.EdgeColor = 'none';
% title('Thickness (m)','FontWeight','normal')
% c1 = [0.7*ones(1,71); 0:0.01:0.7; 0:0.01:0.7  ].' / 0.7;
% c2 = [ 2:-0.01:0; 2:-0.01:0; 2*ones(1,201);].'/2;
% c = [c1; c2];
% xlim([0 1])
% ylim([0 1])
% ylabel('y (km)')
% clim([0.3 3])
% colormap(ax1,c)
% colorbar
% axis square
% set(gca, 'FontSize',12)
% 
% 
% nexttile
% pl2 = pcolor(X/1000+0.45,Y/1000+0.45,abs(mu),'FaceColor','interp');
% pl2.EdgeColor = 'none';
% title('|\mu|','FontWeight','normal')
% colorbar
% xlim([0 1])
% ylim([0 1])
% %clim([-3 3])
% axis square
% set(gca, 'FontSize',12)
% 
% annotation('arrow',[0.905 0.85],[0.11 0.2])
% annotation('arrow',[0.905-0.44 0.85-0.44],[0.11 0.2])
% 
% 
% nexttile
% pl2 = pcolor(X/1000+0.45,Y/1000+0.45,real(phi_tot),'FaceColor','interp');
% pl2.EdgeColor = 'none';
% title('Re(\phi_z)','FontWeight','normal')
% colorbar
% xlim([0 1])
% ylim([0 1])
% ylabel('y (km)')
% xlabel('x (km)')
% %clim([-3 3])
% axis square
% set(gca, 'FontSize',12)
% xlabel('x (km)')
% 
% 
% nexttile
% pl2 = pcolor(X/1000+0.45,Y/1000+0.45,abs(phi_n_tot),'FaceColor','interp');
% pl2.EdgeColor = 'none';
% title('|\phi_z|','FontWeight','normal')
% colorbar
% xlim([0 1])
% ylim([0 1])
% xlabel('x (km)')
% %clim([-3 3])
% axis square
% set(gca, 'FontSize',12)
% xlabel('x (km)')
% 
% fontname(gcf, 'CMU Serif')
% 
% return

%%

saveas(gcf,'spiral15.fig','fig')
exportgraphics(gcf,'spiral15.pdf','ContentType','vector','Resolution',6000)