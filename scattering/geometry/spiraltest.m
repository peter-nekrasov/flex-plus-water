s = 0:0.01:9*pi;
a = 0.0001;
fxt = @(s) -a*s.^3.*cos(-s);
fyt = @(s) -a*s.^3.*sin(-s);
dfxt = @(s) -3*a*s.^2.*cos(-s) - a*s.^3.*sin(-s) ;
dfyt = @(s) -3*a*s.^2.*sin(-s) + a*s.^3.*cos(-s) ;
plot(fxt(s),fyt(s));

rlines = [xs(1),ys(1),xs(end),ys(end)].';
dsig = 0.1;
rlengt = 1;
nlegs = 500; % ceil(4*rlengt/dsig);

quadnodesx = {};
quadnodesy = {};
quadweighs = {};
x = [];
y = [];
w = [];
for ii=1:numel(nlegs)
    [xl,wl,~,~] = lege.exps(nlegs(ii));
    xl = 8*pi*(xl+1)/2 + 3*pi;
    quadweighs{ii} = wl*rlengt(ii)/2.*sqrt(dfxt(xl).^2 + dfyt(xl).^2);
    quadnodesx{ii} = fxt(xl);
    quadnodesy{ii} = fyt(xl);
    x = [x;quadnodesx{ii}];
    y = [y;quadnodesy{ii}];
    w = [w;quadweighs{ii}];
end

xmin = -4;
xmax = 5;
ymin = -4;
ymax = 5;

h = 0.025;
xs = xmin:h:xmax;
ys = ymin:h:ymax;
[XG,YG] = ndgrid(xs,ys);
sz = size(XG);
XG = XG(:);
YG = YG(:);

[xtar,xsrc] = ndgrid(XG,x);
[ytar,ysrc] = ndgrid(YG,y);
amat = exp(-((xtar-xsrc).^2+(ytar-ysrc).^2)/(2*dsig^2)) ;
amatx = amat .* (- (xtar - xsrc) /  dsig^2) ; 
amaty = amat .* (- (ytar - ysrc) /  dsig^2) ; 
amatxx = amat .* (- 1 /  dsig^2) + amatx .* (- (xtar - xsrc) /  dsig^2);
amatxy = amaty .* (- (xtar - xsrc) /  dsig^2) ; 
amatyy = amat .* (- 1 /  dsig^2) + amaty .* (- (ytar - ysrc) /  dsig^2);

hvec = amat*w;
H = reshape(hvec,sz).';
XG = reshape(XG,sz).';
YG = reshape(YG,sz).';


hvec = amat*w;
H = reshape(hvec,sz).';

hxvec = amatx*w;
Hx = reshape(hxvec,sz).';

hyvec = amaty*w;
Hy = reshape(hyvec,sz).';

hxxvec = amatxx*w;
Hxx = reshape(hxxvec,sz).';

hxyvec = amatxy*w;
Hxy = reshape(hxyvec,sz).';

hyyvec = amatyy*w;
Hyy = reshape(hyyvec,sz).';

H0 = 2;
amp = -1;
freq = 1;
gbar = 0;

const = amp/max(H(:));
H = H0+ H*const;
Hx = Hx*const;
Hy = Hy*const;
Hxx = Hxx*const;
Hxy = Hxy*const;
Hyy = Hyy*const;

alpha = E*H.^3/(12*(1-nu^2));
beta = (rhoi*H*freq^2 - rhow*g);

abar = alpha - a0;
bbar = beta - b0;

ax = 3*E*H.^2.*Hx/(12*(1-nu^2));
ay = 3*E*H.^2.*Hy/(12*(1-nu^2));

axx = 6*E*H.*Hx.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hxx/(12*(1-nu^2));
axy = 6*E*H.*Hx.*Hy/(12*(1-nu^2)) + 3*E*H.^2.*Hxy/(12*(1-nu^2));
ayy = 6*E*H.*Hy.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hyy/(12*(1-nu^2));

coefs = {a0,abar,b0,bbar,g0,gbar,ax,ay,axx,axy,ayy,nu}; 
Hs = {H,Hx,Hy,Hxx,Hxy,Hyy};

[X,Y] = meshgrid(xs,ys);

pl1 = pcolor(X,Y,H);
pl1.EdgeColor = 'none';
colorbar


