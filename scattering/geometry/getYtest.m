% checking the derivatives for Y geometry

h = 0.01;

ymin = -0.5;
ymax =  6.5;
xmin = -3.5;
xmax =  3.5;

xs = -xmin:h:xmax;
ys = -ymin:h:ymax;

[X,Y,coefs,Hs] = get_Y(xmin,xmax,ymin,ymax,h,-0.995,0.1);

H = Hs{1};
Hx = Hs{2};
Hy = Hs{3};
Hxx = Hs{4};
Hxy = Hs{5};
Hyy = Hs{6};

ind = find((X == -1) & (Y == 4.5));
ind = ind(1);
[i,j] = ind2sub(size(H),ind);

% stencils

d1 = [1/280; -4/105; 1/5; -4/5; 0; 4/5; -1/5; 4/105; -1/280] / h;

d2 = zeros(9, 1);
d2(1) = -1/560;
d2(2) = 8/315;
d2(3) = -1/5;
d2(4) = 8/5;
d2(5) = -205/72;
d2(6) = 8/5;
d2(7) = -1/5;
d2(8) = 8/315;
d2(9) = -1/560;
d2 = d2 / h^2;

(Hy(i,j) - H(i-4:i+4,j).'*d1) / Hy(i,j)
(Hx(i,j) - H(i,j-4:j+4)*d1) / Hx(i,j)


(Hxx(i,j) - H(i,j-4:j+4)*d2) / Hxx(i,j)
(Hyy(i,j) - H(i-4:i+4,j).'*d2) / Hyy(i,j)


dxy = d1*d1.';

(Hxy(i,j) - sum(dxy.*H(i-4:i+4,j-4:j+4),'all')) / Hxy(i,j)