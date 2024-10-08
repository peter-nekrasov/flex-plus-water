r = 0:0.01:60;
rhoj = 0.304723573570608 + 1.13327427908085i;
zt = r*rhoj;
x1 = besselh(0, zt);

src = [0; 0];
targ = [zt / rhoj; 0*r];
[h0, ~] = helmdiffgreen(rhoj,src,targ);
h0 = -4i*h0.' + 4i/(2*pi)*log(zt / rhoj);

tiledlayout(1,3);
nexttile
plot(r, real(x1),r,imag(x1))
legend('real', 'imag')
title('besselh()')

nexttile
plot(r, real(h0),r,imag(h0))
legend('real', 'imag')
title('helmdiffgreen.m')

nexttile
plot(r, real(abs(x1-h0)))