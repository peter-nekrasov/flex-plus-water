
rlines = [0,0,0,4-0.2*cos(pi/4);0.1,4,2,6;-0.1,4,-2,6].';
dsig = 0.2;
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

ymin = -0.5;
ymax =  6.5;
xmin = -3.5;
xmax =  3.5;

xs = xmin:0.01:xmax;
ys = ymin:0.01:ymax;
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
H = reshape(hvec,sz);

hxvec = amatx*w;
Hx = reshape(hxvec,sz);

hyvec = amaty*w;
Hy = reshape(hyvec,sz);

hxxvec = amatxx*w;
Hxx = reshape(hxxvec,sz);

hxyvec = amatxy*w;
Hxy = reshape(hxyvec,sz);

hyyvec = amatyy*w;
Hyy = reshape(hyyvec,sz);

[X,Y] = meshgrid(xs,ys);

figure(1);
tiledlayout(2,3);
nexttile
p1 = pcolor(X,Y,H);
p1.EdgeColor = 'none';
colorbar

nexttile
p1 = pcolor(X,Y,Hx);
p1.EdgeColor = 'none';
colorbar

nexttile
p1 = pcolor(X,Y,Hy);
p1.EdgeColor = 'none';
colorbar

nexttile
p1 = pcolor(X,Y,Hxx);
p1.EdgeColor = 'none';
colorbar

nexttile
p1 = pcolor(X,Y,Hxy);
p1.EdgeColor = 'none';
colorbar

nexttile
p1 = pcolor(X,Y,Hyy);
p1.EdgeColor = 'none';
colorbar