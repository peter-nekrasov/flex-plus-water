
rlines = [0,0,0,4;0,4,2,6;0,4,-2,6].';
dsig = 0.1;
rlengt = sqrt((rlines(3,:)-rlines(1,:)).^2 + (rlines(4,:)-rlines(2,:)).^2);
nlegs = ceil(4*rlengt/dsig);
﻿
quadnodesx = {};
quadnodesy = {};
quadweighs = {};
x = [];
y = [];
w = [];
for ii=1:numel(nlegs)
    [xl,wl,~,~] = lege.exps(nlegs(ii));
    quadweighs{ii} = wl*rlengt(ii)/2;
    quadnodesx{ii} = (rlines(3,ii)-rlines(1,ii))*(xl+1)/2+rlines(1,ii);
    quadnodesy{ii} = (rlines(4,ii)-rlines(2,ii))*(xl+1)/2+rlines(2,ii);
    x = [x;quadnodesx{ii}];
    y = [y;quadnodesy{ii}];
    w = [w;quadweighs{ii}];
end
﻿
﻿
ymin = -0.5;
ymax =  6.5;
xmin = -3.5;
xmax =  3.5;
﻿
xs = xmin:0.01:xmax;
ys = ymin:0.01:ymax;
[XG,YG] = ndgrid(xs,ys);
sz = size(XG);
XG = XG(:);
YG = YG(:);
﻿
[xtar,xsrc] = ndgrid(XG,x);
[ytar,ysrc] = ndgrid(YG,y);
amat = exp(-((xtar-xsrc).^2+(ytar-ysrc).^2)/(2*dsig^2));
hvec = amat*w;
﻿
H = reshape(hvec,sz);