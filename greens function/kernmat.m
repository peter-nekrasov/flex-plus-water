function [intmat] = kernmat(intfun,h,xmax,nt)
%KERNMAT forms the square matrix for a kernel using the 
%       center-corrected trapezoid rule outlined in Duan & Rokhlin (2008).
%
% input:
%
% intfun - @(x,y) function handle of integral operator K(x,y)
% h - float, grid spacing 
% xmax - float, box size. kernel will be evaluated for x,y in [-xmax,xmax]
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


if nargin < 4
    nt = 5;
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

[xs, ys] = meshgrid(-xmax:h:xmax);
[gmat,~] = intfun(xs,ys);

[xs, ys] = meshgrid(-(lev-1):(lev-1));
rs = sqrt(xs.^2 + ys.^2);
gmat = ;

pts = sortrows([rs(:) xs(:) ys(:)]);
pts = pts(1:nt,2:3).*h;

pts = pts((pts(:,1) >= 0) & (pts(:,2) >= 0),:);
xpts = pts(:,1);
ypts = pts(:,2);
[npt,~] = size(pts);

isevenx = ( fun([1,0]) == fun(-1,0) );
iseveny = ( fun([0,-1]) == fun(0,-1) );

if isevenx
    xpow = 0;
else
    xpow = 1;
end

if iseveny
    ypow = 0;
else 
    ypow = 1;
end

xyswitch = 1;

for pow = 1:npt
    lhs(pow,:) = (xpts.^xpow).*(ypts.^ypow);

%    rhs = integral2(intfun)*exp() - mat;


    if (xyswitch < 0)
        xpow = xpow+1;
    else
        ypow = ypow+1;
    end

    xyswitch = -xyswitch;
end