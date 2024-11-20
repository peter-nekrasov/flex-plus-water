% shallow waves

ks = -0.1:0.001:0.1;
h = 1;
g = 9.8;
ws = sqrt(g*ks.*tanh(ks*h));

plot(ws,ks,'b',-ws,ks,'b')
hold on

% deep waves

ws = sqrt(g*abs(ks));

plot(ws,ks,'g',-ws,ks,'g')
hold on


% hydroelastic deep waves


alpha = 10^11;
rhow = 1025;
rhoi = 950;

ws = sqrt(alpha*abs(ks).^5 + rhow*g*abs(ks))./(rhow + rhoi*abs(ks));

plot(ws,ks,'k',-ws,ks,'k')
hold on

xlabel('\omega');
ylabel('k');

legend('shallow waves','deep waves','flexural waves (deep)')
