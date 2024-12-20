function [X,Y,coefs,Hs] = get_Y(xmin,xmax,ymin,ymax,h,amp,freq)

E = 7*10^9;
nu = 0.33;
H0 = 1;
rhow = 1025;
rhoi = 917;
g = 9.8;

a0 = E*H0^3/(12*(1-nu^2));
b0 = (rhoi*H0*freq^2 - rhow*g);
g0 = -freq^2*rhow;
gbar = 0;

rlines = [0,0,0,40-2*cos(pi/4);1,40,20,60;-1,40,-20,60].';
dsig = 2;
rlengt = sqrt((rlines(3,:)-rlines(1,:)).^2 + (rlines(4,:)-rlines(2,:)).^2);
nlegs = ceil(4*rlengt/dsig);

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

end