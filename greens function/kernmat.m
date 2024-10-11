function intmat = kernmat(intfun,h,ymax,nt)
%KERNMAT forms the square matrix for a kernel using the 
%       center-corrected trapezoid rule outlined in Duan & Rokhlin (2008).
%
% input:
%
% intfun - @(y1,y2) function handle of integral operator K(x1,x2,y1,y2) 
%               where (x1,x2) = (0,0).
% h - float, grid spacing 
% ymax - float, box size with y1, y2 in [-ymax, ymax].
%
% optional input: 
%
% nt - int, number of points in the stencil. default: nt = 5. nt must be a
%   centered square number. 
%           - i.e. nt = 1, fourth order
%           - nt = 5, sixth order
%           - nt = 13, eighth order
%           - nt = 25, tenth order, etc. 
%
% output: 
%
% - intmat is the center-corrected integral operator 


[X, Y] = meshgrid(-ymax:h:ymax);
[gmat] = intfun(X,Y);
gmat(isnan(gmat)) = 0;

if nargin < 4
    nt = 5;
elseif nt == 0
    intmat = gmat*h^2;
    return
else
    lev = 0;
    numb = 0;

    while numb < nt
        numb = lev^2 + (lev+1)^2;
        lev = lev + 1;
    end

    if (nt < 1) || (numb ~= nt)
        error('nt must be a centered square number')
    end
end

lhs = zeros(nt);
rhs = zeros(nt,1);

[xs, ys] = meshgrid(-(lev-1):(lev-1));
rs = sqrt(xs.^2 + ys.^2);

pts = sortrows([rs(:) xs(:) ys(:)]);
iis = pts(1:nt,2);
jjs = pts(1:nt,3);
xpts = pts(1:nt,2)*h;
ypts = pts(1:nt,3)*h;


if intfun(1,1) == -intfun(-1,1)
    xodd = 1;
elseif intfun(1,1) == intfun(-1,1)
    xodd = 0;
else
    xodd = 0.5;
end

if intfun(1,1) == -intfun(1,-1)
    yodd = 1;
elseif intfun(1,1) == intfun(1,-1)
    yodd = 0;
else
    yodd = 0.5;
end

i = 1;
z = 0;
while i <= nt
    [xpow, ypow] = cantor(z);
    row = (xpts.^xpow).*(ypts.^ypow).*exp(-1/(25*h^2)*(xpts.^2+ypts.^2));
    lhs(i,:) = row.';
    if rank(lhs) == i
        
        if (mod(xodd + xpow,2) == 1) || (mod(yodd + ypow,2) == 1) 
            rhs(i) = 0;
        else
            ftest = exp(-1/(25*h^2)*(X.^2+Y.^2)).*(X.^xpow).*(Y.^ypow);
            integrand = @(x,y) intfun(x,y).*exp(-1/(25*h^2)*(x.^2+y.^2)).*(x.^xpow).*(y.^ypow);
            ymax = min(30*h,ymax);
            rhs(i) = integral2(integrand,-ymax,ymax,-ymax,ymax,"RelTol",1e-13,'AbsTol',0)/h^2 - sum(gmat.*ftest,'all');
        end
        
        i = i + 1;
    end
    z = z + 1;
end

corrections = lhs \ rhs;

ind = intersect(find(X == 0), find(Y == 0));
[centx, centy] = ind2sub(size(X),ind);


for i = 1:nt
    gmat(iis(i)+centx,jjs(i)+centy) = gmat(iis(i)+centx,jjs(i)+centy) + corrections(i);
end

intmat = gmat*h^2;

end


function [x,y] = cantor(z)
%
% Gets Cantor pairing from a positive integer z
%
    w = floor((sqrt(8*z+1)-1)/2);
    t = (w^2+w)/2;
    y = z - t;
    x = w-y;
end