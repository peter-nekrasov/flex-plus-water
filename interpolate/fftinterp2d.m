% 2d interpolation using fft

% generate some random data

l = 200;
h = 10;
xs = (-l/2:h:l/2-h);
N = length(xs);
[X,Y] = meshgrid(xs);

Z = exp(2i*pi*X/l) + 0.5*exp(6i*pi*Y/l) + 0.25*exp(4i*pi*X/l).*exp(4i*pi*Y/l);

figure(1);
tiledlayout(2,2)

nexttile
s = pcolor(X,Y,real(Z));
s.EdgeColor = 'none';
title('data - real part')
colorbar

nexttile 
s = pcolor(X,Y,imag(Z));
s.EdgeColor = 'none';
title('data - imaginary part')
colorbar

% perform fft

z = (fft2(Z));

% interpolation - zero padding

hp = 1;
px = (-l/2:hp:l/2-hp);
Np = length(px);
[Xp,Yp] = meshgrid(px);
Zp = ifft2(z,Np,Np)*Np^2/N^2;

% plot

nexttile
s = pcolor(Xp,Yp,real(Zp));
s.EdgeColor = 'none';
title('interpolation - real part')
colorbar

nexttile 
s = pcolor(Xp,Yp,imag(Zp));
s.EdgeColor = 'none';
title('interpolation - imaginary part')
colorbar

err = max(abs(Zp(ismember(Xp,X) & ismember(Yp,Y)) - Z(:)))
