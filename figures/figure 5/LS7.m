%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for plane wave 
% scattering of flexural-gravity waves with rolls
%
%%%%%

h = 20;

x1 = -1E3;
x2 = 5E3;

y1 = -3E3;
y2 = 3E3;

xs = x1:h:x2;
ys = y1:h:y2;
xl = 2*x1:h:2*x2;
yl = 2*y1:h:2*y2;
[~,n] = size(xs);
[X,Y] = meshgrid(xs,ys);
[XL,YL] = meshgrid(xl,yl);

srcinfo = []; srcinfo.r = [x1;0];
targinfo = []; targinfo.r = [X(:) Y(:)].';

freqs = 0.0025:0.0025:4;
ks = freqs*0;
Rs = freqs*0;
Ts = freqs*0;

mus = zeros(length(freqs), length(X(:)));
phis = zeros(length(freqs), length(X(:)));
phizs = zeros(length(freqs), length(X(:)));
flags = zeros(length(freqs));
relreses = zeros(length(freqs));
iters = zeros(length(freqs),2);
resvecs = cell(length(freqs));

mu_pre = zeros(length(X(:)),1);

for ii = 1:length(freqs)

    w = freqs(ii);

    fprintf('Frequency: %.4f ', w);

    % if w > 0.5
    %     h = 12.5;
    % end
    
    [coefs, H] = rolls(X,Y,0,40E2,-25E2,25E2,0.75,333.3,w); % remove gbar from coefs vector

    a0 = coefs{1}; 
    b0 = coefs{3}; 
    g0 = coefs{5}; 
    
    % Finding positive real roots
    [rts,ejs] = find_roots(b0 / a0, g0 / a0);
    k = rts((imag(rts) == 0) & (real(rts) > 0));
    ks(ii) = k;
    ejs = ejs/a0;
    
    % RHS (Incident field)
    kerns = green(srcinfo.r,targinfo.r,rts,ejs);
    phizinc = reshape(kerns{1},size(X));
    phiinc = reshape(kerns{4},size(X));
    [rhs_vec] = get_rhs_vec2(coefs,kerns);
    
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
    [mu, flag, relres, iter, resvec] = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,30,1e-6,500,[],[],mu_pre);

    mu_pre = mu;
    if flag 
        fprintf('GMRES failed with relative residual %d \n', relres);
    else 
        fprintf('GMRES succeeded. \n');    
    end

    mu = reshape(mu, size(X));
    
    [phi, phi_z] = sol_eval_fft(mu,evalkerns);
    phi = reshape(phi,size(X));
    phi_z = reshape(phi_z,size(X));

    phi_tot = phi + phiinc;
    phi_z_tot = phi_z + phizinc;

    mus(ii,:) = mu(:);
    phis(ii,:) = phi_tot(:);
    phizs(ii,:) = phi_z_tot(:);
    flags(ii) = flag;
    relreses(ii) = relres;
    resvecs{ii} = resvec;

    f = figure(1);
    set(f,'Position',[1 1200 1332 300])
    t = tiledlayout(1,3);
    nexttile
    pc = pcolor(X,Y,real(phi_tot));
    pc.EdgeColor = 'none';
    title('Re(\phi)')
    colorbar

    nexttile
    pc = pcolor(X,Y,real(phi_z_tot));
    pc.EdgeColor = 'none';
    title('Re(\partial_z\phi)')
    colorbar

    nexttile
    pc = pcolor(X,Y,real(mu));
    pc.EdgeColor = 'none';
    title('Re(\mu)')
    colorbar
    drawnow
    
    Tind = find((X == 4500) & (Y == 0));
    Rind = find((X == (-500)) & (Y == 0));

    Ts(ii) = abs(phi_tot(Tind));
    Rs(ii) = abs(phi(Rind));

end

%%

figure(2); clf
plot(freqs(1:ii-1),imag(phis(1:ii-1,2000)))


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
clim([0 max(abs(phi_z_tot(:)))*0.75])
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
pc = pcolor(X,Y,abs(phi_z_tot));
clim([0 max(abs(phi_z_tot(:)))*0.75])
pc.EdgeColor = 'none';
%title('|\phi_n|')
axis off



%%

figure(1);

tiledlayout(2,1,"TileSpacing","tight")

% nexttile
% C = imread('/Users/peter/Downloads/Ward_Hunt_Island,_Ice_Shelf_02.jpg');
% image(C)
% axis off


h = 10;

x1 = -1E3;
x2 = 5E3;

y1 = -3E3;
y2 = 3E3; 

xs = x1:h:x2;
ys = y1:h:y2;
xl = 2*x1:h:2*x2;
yl = 2*y1:h:2*y2;
[~,n] = size(xs);
[X,Y] = meshgrid(xs,ys);
[XL,YL] = meshgrid(xl,yl);
[coefs, H] = rolls(X,Y,0,40E2,-25E2,25E2,0.75,333.3,1); % remove gbar from coefs vector


c = [0.95:-0.01:0.3 ; 0.95:-0.01:0.3; 0.95:-0.01:0.3 ].';


nexttile
s = pcolor(X/1000,Y/1000,H);
s.EdgeColor = 'none';
%xlim([0.4 5])
%ylim([-2 2])
title('Thickness (m)','FontWeight','normal')
colormap(c);
colorbar
hold on
%plot((-3:0.01:3)*0,-3:0.01:3,'k--','LineWidth',1)
%hold on
scatter(-0.5,0,50,'filled')
hold on
scatter(4.5,0,50,'filled')
xlabel('x (km)')
ylabel('y (km)')
annotation('arrow',[0.14 0.18],[0.69 0.69])
annotation('arrow',[0.14 0.18],[0.74 0.74])
annotation('arrow',[0.14 0.18],[0.64 0.64])
annotation('arrow',[0.14 0.18],[0.79 0.79])
annotation('arrow',[0.14 0.18],[0.84 0.84])


%quiver([-0.9 -0.9 -0.9 -0.9 -0.9],[-2 -1 0 1 2],[0.2 0.2 0.2 0.2 0.2],[0 0 0 0 0],'off',"Color","black","ShowArrowHead","on",LineWidth=3)

set(gca, 'FontSize',11)


Rs1 = load("Rs1.mat").Rs(1:651);
Rs2 = load("Rs2.mat").Rs(1:25);
Rs3 = load("Rs3.mat").Rs;
Rs4 = load("Rs4.mat").Rs;

Ts1 = load("Ts1.mat").Ts(1:651);
Ts2 = load("Ts2.mat").Ts(1:25);
Ts3 = load("Ts3.mat").Ts;
Ts4 = load("Ts4.mat").Ts;

freqs1 = load("freqs1.mat").freqs(1:651);
freqs2 = load("freqs2.mat").freqs(1:25);
freqs3 = load("freqs3.mat").freqs;
freqs4 = load("freqs4.mat").freqs;

ks1 = load("ks1.mat").ks(1:651);
ks2 = load("ks2.mat").ks(1:25);
ks3 = load("ks3.mat").ks;
ks4 = load("ks4.mat").ks;

Rs = [Rs1 Rs2 Rs3 Rs4];
Ts = [Ts1 Ts2 Ts3 Ts4];
freqs = [freqs1 freqs2 freqs3 freqs4];
ks = [ks1 ks2 ks3 ks4];

nexttile
p = plot(ks,Rs,ks,Ts,'LineWidth',1)

%legend('Reflected','Transmitted','Location','best')
xlabel('k (m^{-1})')
xlim([min(ks) max(ks)])
ylabel('|\phi|')
hold on

plot((0:0.01:1.2)*0+0.0235,0:0.01:1.2,'k--','LineWidth',0.7)
hold on

plot((0:0.01:1.2)*0+0.0265,0:0.01:1.2,'k--','LineWidth',0.7)

legend([p(1) p(2)],{'Reflected','Transmitted'},'Location','east')

set(gca, 'FontSize',11)

fontname(gcf, 'CMU Serif')



% nexttile
% pc = pcolor(X/1000,Y/1000,abs(phi_n_tot));
% ylim([-6 6])
% pc.EdgeColor = 'none';
% %title('|\phi_n|')
% colorbar
% colormap(gca,"default")
% drawnow
% set(gca, 'FontSize',12)


fontname(gcf, 'CMU Serif')

%% 
saveas(gcf,'rollfig.fig','fig')
exportgraphics(gcf,'rollfig.pdf','ContentType','vector')

%%

f= figure(1);
tiledlayout(1,2,'TileSpacing','tight')

X = load('X1.mat').X;
Y = load('Y1.mat').Y;
phi1 = load('phi_n_tot2.mat').phi_n_tot;
phi2 = load('phi_n_tot1.mat').phi_n_tot;

nexttile
pc = pcolor(X/1000,Y/1000,abs(phi1));
ylim([-6 6])
pc.EdgeColor = 'none';
title('k = 0.0235','FontWeight','normal')
clim([0 max(abs(phi2(:)))])
drawnow
set(gca, 'FontSize',12)
xlabel('x (km)')
ylabel('y (km)')
axis square

nexttile
pc = pcolor(X/1000,Y/1000,abs(phi2));
ylim([-6 6])
pc.EdgeColor = 'none';
title('k = 0.0265','FontWeight','normal')
colorbar
clim([0 max(abs(phi2(:)))])
drawnow
set(gca, 'FontSize',12)
xlabel('x (km)')
fontname(gcf, 'CMU Serif')
axis square

%% 
saveas(gcf,'rollfig2.fig','fig')
exportgraphics(gcf,'rollfig3.pdf','ContentType','vector')