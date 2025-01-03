
function [proxy,pw] = proxy_ann_pts2(p0,p,q,width)

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

theta = (0:(p0-1))*2*pi/p0;
T = theta(:);
pw0 = ones(p0,1)*(2.0*pi)/p0;
P = pw0(:);
proxy = [proxy,[1.5*cos(T).';1.5*sin(T).']];
pw = [pw,pw0.'];

theta = (0:(p0-1))*2*pi/p0;
T = theta(:);
pw0 = ones(p0,1)*(2.0*pi)/p0;
P = pw0(:);
proxy = [proxy,[(1.5+width)*cos(T).';(1.5+width)*sin(T).']];
pw = [pw,pw0.'];

end