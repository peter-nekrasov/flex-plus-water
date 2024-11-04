% test convolutions

addpath(genpath('..'))
h = 0.125;

xs = -15:h:15;
xl = -2*15:h:2*15;
[~,n] = size(xs);
[X,Y] = meshgrid(xs);
[XL,YL] = meshgrid(xl);

targ = [2; 2];
src = [X(:).'; Y(:).'];
a0 = 1;
gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);




dens = @(x,y) x.*exp(-(x.^2+y.^2)/(10));
truev =  -0.089213311677999 + 0.116604914243712i;% integral2(@(x,y) dens(x,y).*greenfac(x,y),-20,20,-20,20,"AbsTol",0,"RelTol",10E-16);
[inds,corrs] = get_correct(rts,ejs,h,a0);
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),inds,corrs);
kern = kerns{3};
kern = kern(:,:,1);
d1 = dens(X,Y);
abs(truev - sum(kern(:).*d1(:)))


%%

ind = find((XL == 0) & (YL ==0));
sz = size(XL);

src = [0; 0];
targ = [XL(:).'; YL(:).'];
kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),inds,corrs);
kerns = gen_fft_kerns(kerns,sz,ind);
kernhat = kerns{3};
kernhat = kernhat(:,:,1);

d1aug = [d1, zeros(n,n-1); zeros(n-1,n), zeros(n-1)];
d1hat = fft2(d1aug);

prodmat = ifft2(kernhat.*d1hat);
outmat = prodmat(1:n, 1:n);

ind2 = find((X == 2) & (Y == 2));
abs(outmat(ind2) - truev)
