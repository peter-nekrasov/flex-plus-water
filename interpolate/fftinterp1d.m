% 1d interpolation using fft

% generate some random data

l = 200;
h1 = 20;
xs = -l/2:h1:l/2;
N = length(xs);

func = @(x) exp(-(x).^2/100);
ys = func(xs); % exp(2i*pi*xs/l) + 0.5*exp(6i*pi*xs/l);

plot(xs,real(ys),'o','DisplayName','samples')
hold on

% plot(xs,imag(ys),'bo')
% 
% ys = interpft(ys,1000);
% xs = (0:1000-1)/1000*l;
% 
% plot(xs,real(ys),'-','DisplayName','samples')
% xlim([0 l])
% hold on


% perform fft

z = fftshift(fft(ys));

% interpolation method 1 - summation over fourier 

% Np = 80;
% px = (0:Np-1)/Np*l;
% k = -floor(N/2):floor(N/2);
% py = z/N*exp(2i*pi*k'*px/l);

% plot(px,real(py),'r-','DisplayName','real')



% interpolation method 2 - zero padding

h = 0.25;
px = -l/2:h:l/2+h1-h;
Np = length(px);
npad = Np-N;
z = [zeros(1,ceil(npad/2)) z zeros(1,floor(npad/2))];
z = ifftshift(z);
py = ifft(z)*Np/N;

% plot

plot(px,real(py),'DisplayName','FFT interpolation')
hold on

plot(px,func(px),'DisplayName','original function')
% plot(px,imag(py),'b-','DisplayName','imaginary')

legend

err = max(abs(py(ismember(px,xs)) - ys))
