mu5 = load('mu_Y_005.mat').sol2;
%mu2 = load('mu_Y_002.mat').sol2;
%mu1 = load('mu_Y_001.mat').sol2;
mu05 = load('mu_Y_0005.mat').sol2;

phi1 = load('phi_Y_001.mat').phi;
phi2 = load('phi_Y_002.mat').phi;
phi5 = load('phi_Y_005.mat').phi;
phi05 = load('phi_Y_0005.mat').phi;

phin1 = load('phin_Y_001.mat').phi_n;
phin2 = load('phin_Y_002.mat').phi_n;
phin5 = load('phin_Y_005.mat').phi_n;
phin05 = load('phin_Y_0005.mat').phi_n;


X = load('X_Y.mat').X;
Y = load('Y_Y.mat').Y;

inind = find((X == -8) & (Y == 0.5));
outind = find((X == -8) & (Y == 20));
%bind = find((X == 100) & (Y == 0));

hs = [0.005 0.01 0.02 0.05];
ins = [ phi05(inind) phi1(inind) phi2(inind) phi5(inind)];
outs = [ phi05(outind) phi1(outind) phi2(outind) phi5(outind)];
% bs = [ phi05(bind) phi1(bind) phi2(bind) phi5(bind)];

figure(1);
tiledlayout(1,1);
nexttile
plot(hs,real(ins),'o-')
title('inside solution')
xlabel('min thickness')
ylabel('real(\phi)')

% nexttile
% plot(hs,abs(outs),'o-')
% title('outside solution')
% xlabel('min thickness')

% nexttile
% plot(hs,bs,'o-')
% xlim([0 0.01])
% title('boundary solution')
% xlabel('min thickness')

%% finding y intercept

A = [1 0.001; 1 0.002];
rhs = [phi1(:).'; phi2(:).'];

x = A \ rhs;

phi0 = x(1,:);
phi0 = reshape(phi0,size(X));

rhs = [phin1(:).'; phin2(:).'];

x = A \ rhs;

phin0 = x(1,:);
phin0 = reshape(phin0,size(X));

rhs = [mu1(:).'; mu2(:).'];

x = A \ rhs;

mu0 = x(1,:);
mu0 = reshape(mu0,size(X));

k1 = 0.059374132319434; k2 = 0;
phiinc = exp(1i*k1*X+1i*k2*Y);

phi0 = phi0 + phiinc;
phin0 = phin0 + k*phiinc;

figure(2);
tt = tiledlayout(1,3);
nexttile
p1 = pcolor(X,Y,real(phin0));
p1.EdgeColor = 'none';
title('re(\phi_n)')
colorbar
axis square

nexttile
p1 = pcolor(X,Y,abs(phin0));
p1.EdgeColor = 'none';
title('|\phi_n|')
colorbar
axis square

nexttile
p1 = pcolor(X,Y,abs(mu0));
p1.EdgeColor = 'none';
title('|\mu|')
colorbar
axis square

title(tt,'Richardson extrapolation')