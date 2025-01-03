
function [proxy,pw] = proxy_ann_pts(p,q,width)

[gq,wq,~,~] = lege.exps(q);
gq = (gq+1)*width/2 + 1.5;
wq = wq*width/2;

theta = (0:(p-1))*2*pi/p;
pw = ones(p,1)*(2.0*pi)/p;
rads  = gq;

[R,T] = meshgrid(rads,theta);
[WR,WT] = meshgrid(wq,pw);
R = R(:);
T = T(:);
W = WR(:).*WT(:).*R;

proxy = [(R.*cos(T)).';(R.*sin(T)).'];
pw = W.';

end