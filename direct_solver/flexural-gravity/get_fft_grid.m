function [src,targ,ind,sz] = get_fft_grid(N,L)

xl = 2*L*(-(N-1):(N-1))/(N-1);
[XL,YL] = meshgrid(xl);

src = [0;0];
targ = [XL(:).'; YL(:).'];

ind = find((XL == 0) & (YL ==0));
sz = size(XL);

end