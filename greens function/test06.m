h = 0.2;
xs = -15:h:15;
[X,Y] = meshgrid(xs);
targ = [X(:).'; Y(:).'];
src = [0; 0];
gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);

akern = green(src,targ,rts,ejs);
bkern = green(targ,src,rts,ejs);

a1 = akern{1};
a2 = akern{2};
a3 = akern{3};
a4 = akern{4};

b1 = bkern{1};
b2 = bkern{2};
b3 = bkern{3};
b4 = bkern{4};

err = max(abs(a1 - b1.'))
err2 = max(abs(a2(:,:,1) - b2(:,:,1).'))
err3 = max(abs(a2(:,:,2) - b2(:,:,2).'))
err4 = max(abs(a2(:,:,3) - b2(:,:,3).'))
err5 = max(abs(a3(:,:,1) + b3(:,:,1).'))
err6 = max(abs(a3(:,:,2) + b3(:,:,2).'))
err7 = max(abs(a3(:,:,3) + b3(:,:,3).'))
err8 = max(abs(a3(:,:,4) + b3(:,:,4).'))
