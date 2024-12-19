mu1 = load('mu0001.mat').mu;
mu2 = load('mu0002.mat').mu;
mu4 = load('mu0004.mat').mu;
mu8 = load('mu0008.mat').mu;

phi1 = load('phi0001.mat').phi;
phi2 = load('phi0002.mat').phi;
phi4 = load('phi0004.mat').phi;
phi8 = load('phi0008.mat').phi;

phin1 = load('phin0001.mat').phi_n;
phin2 = load('phin0002.mat').phi_n;
phin4 = load('phin0004.mat').phi_n;
phin8 = load('phin0008.mat').phi_n;

X = load('X.mat').X;
Y = load('Y.mat').Y;

inind = find((X == 0) & (Y == 0));
outind = find((X == 200) & (Y == 200));
bind = find((X == 100) & (Y == 0));

hs = [0.001 0.002 0.004 0.008];
ins = [phi1(inind) phi2(inind) phi4(inind) phi8(inind)];
outs = [phi1(outind) phi2(outind) phi4(outind) phi8(outind)];
bs = [phi1(bind) phi2(bind) phi4(bind) phi8(bind)];

figure(1);
tiledlayout(1,3);
nexttile
plot(hs,ins,'o-')
xlim([0 0.01])
title('inside solution')
xlabel('min thickness')

nexttile
plot(hs,outs,'o-')
xlim([0 0.01])
title('outside solution')
xlabel('min thickness')

nexttile
plot(hs,bs,'o-')
xlim([0 0.01])
title('boundary solution')
xlabel('min thickness')

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