h = 0.2;
xs = -15:h:15;
[X,Y] = meshgrid(xs);
targ = [X(:).'; Y(:).'];
src = [0; 0];

[a1,~,a2,a3] = helmdiffgreen(6,src,targ);
[b1,~,b2,b3] = helmdiffgreen(6,targ,src);

err = max(abs(a1 - b1.'))
err2 = max(abs(a2(:,:,1) - b2(:,:,1).'))
err3 = max(abs(a2(:,:,2) - b2(:,:,2).'))
err4 = max(abs(a2(:,:,3) - b2(:,:,3).'))
err5 = max(abs(a3(:,:,1) + b3(:,:,1).'))
err6 = max(abs(a3(:,:,2) + b3(:,:,2).'))
err7 = max(abs(a3(:,:,3) + b3(:,:,3).'))
err8 = max(abs(a3(:,:,4) + b3(:,:,4).'))
