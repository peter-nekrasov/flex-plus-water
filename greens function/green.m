function [val,grad,hess,gradlap] = green(x,y,beta,gamma,opt)
%
% computes the green's function centered at (x,y) = 0 for the 
% integro-differential equation determined by the roots of the polynomial:
%             z^5 - beta*z + gamma = 0
%
% 
% output : (note the convention is not the same as helmdiffgreen.m)
% - val has the value of the Green's function centered at zero and
%   evaluated at (x,y)
% - grad(:,:,1) has G_{x}, grad(:,:,2) has G_{y}
% - hess(:,:,1) has G_{xx}, hess(:,:,2) has G_{xy}, 
%   hess(:,:,3) has G_{yy}
% - gradlap is the gradient of the Laplacian, namely 
%   gradlap(:,:,1) has G_{xxx} + G_{xyy}, gradlap(:,:,2) has G_{yxx} + G_{yyy}
% 
%
% input:
%
% x - x-coordinates array
% y - y-coordinates array
% beta - coefficient beta in the equation
% gamma - coefficient gamma in the equation
%
% optional input:
%
% opt - string, default: 'green'. 
%         Possible options are:
%         opt = 'green' => Green's function (and derivatives)
%         opt = 'phi' => kernel used to evaluate phi on surface
%

if nargin < 5
    opt = 'green';
end

r = sqrt(x.^2 + y.^2);
sz = size(r);
r = abs(r);
r = r(:).';
x = x(:).';
y = y(:).';

[rts2, ejs] = find_roots(beta,gamma);

if strcmpi(opt,'phi')
    ejs = ejs./rts2;
end

val = 0;
gradx = 0;
grady = 0;
hessxx = 0;
hessxy = 0;
hessyy = 0;
gradlapx = 0;
gradlapy = 0;

src = [0; 0];
targ = [x; y];

for i = 1:5
    
    rhoj = rts2(i);
    ej = ejs(i);

    if angle(rhoj) == 0

       [sk0,gradsk0,hesssk0,gradlapsk0] = struveKdiffgreen(rhoj,src,targ);
       [h0,gradh0,hessh0,thirdh0] = helmdiffgreen(rhoj,src,targ);

       h0(r == 0) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));

       h0 = -4i*h0;
       gradh0 = -4i*gradh0;
       
       h0x = gradh0(:,:,1);
       h0y = gradh0(:,:,2);
       
       h0x(r == 0) = 0;
       h0y(r == 0) = 0;
       
       h0xx = hessh0(:,:,1);
       h0xy = hessh0(:,:,2);
       h0yy = hessh0(:,:,3);

       h0xx(r == 0) = rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-1.5-log(2));
       h0xy(r == 0) = 0;
       h0yy(r == 0) = rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-1.5-log(2));
       
       h0xx = -4i*h0xx;
       h0xy = -4i*h0xy;
       h0yy = -4i*h0yy;

       h0xxx = thirdh0(:,:,1);
       h0yxx = thirdh0(:,:,2);
       h0xyy = thirdh0(:,:,3);
       h0yyy = thirdh0(:,:,4);

       h0xxx(r == 0) = 0;
       h0yxx(r == 0) = 0;
       h0xyy(r == 0) = 0;
       h0yyy(r == 0) = 0;

       h0xxx = -4i*h0xxx;
       h0yxx = -4i*h0yxx;
       h0xyy = -4i*h0xyy;
       h0yyy = -4i*h0yyy;
       
       sk0x = gradsk0(:,:,1);
       sk0y = gradsk0(:,:,2);

       sk0xx = hesssk0(:,:,1);
       sk0xy = hesssk0(:,:,2);
       sk0yy = hesssk0(:,:,3);

       sk0lapx = gradlapsk0(:,:,1);
       sk0lapy = gradlapsk0(:,:,2);

       val = val + ej*rhoj^2*(-sk0 + 2i*h0);

       gradx = gradx + ej*rhoj^2*(-sk0x + 2i*h0x);
       grady = grady + ej*rhoj^2*(-sk0y + 2i*h0y);
       
       hessxx = hessxx + ej*rhoj^2*(-sk0xx + 2i*h0xx);
       hessxy = hessxy + ej*rhoj^2*(-sk0xy + 2i*h0xy);
       hessyy = hessyy + ej*rhoj^2*(-sk0yy + 2i*h0yy);

       gradlapx = gradlapx + ej*rhoj^2*(-sk0lapx + 2i*(h0xxx+h0xyy));
       gradlapy = gradlapy + ej*rhoj^2*(-sk0lapy + 2i*(h0yxx+h0yyy));


    elseif rhoj ~= 0

       [sk0,gradsk0,hesssk0,gradlapsk0] = struveKdiffgreen(-rhoj,src,targ);

       sk0x = gradsk0(:,:,1);
       sk0y = gradsk0(:,:,2);

       sk0xx = hesssk0(:,:,1);
       sk0xy = hesssk0(:,:,2);
       sk0yy = hesssk0(:,:,3);

       sk0lapx = gradlapsk0(:,:,1);
       sk0lapy = gradlapsk0(:,:,2);

       val = val + ej*rhoj^2*sk0;

       gradx = gradx + ej*rhoj^2*sk0x;
       grady = grady + ej*rhoj^2*sk0y;

       hessxx = hessxx + ej*rhoj^2*sk0xx;
       hessxy = hessxy + ej*rhoj^2*sk0xy;
       hessyy = hessyy + ej*rhoj^2*sk0yy;

       gradlapx = gradlapx + ej*rhoj^2*sk0lapx;
       gradlapy = gradlapy + ej*rhoj^2*sk0lapy;

    end

end

val = pi/2*val;
gradx = pi/2*gradx;
grady = pi/2*grady;
hessxx = pi/2*hessxx;
hessxy = pi/2*hessxy;
hessyy = pi/2*hessyy;
gradlapx = pi/2*gradlapx;
gradlapy = pi/2*gradlapy;

val = reshape(val,sz);
gradx = reshape(gradx,sz);
grady = reshape(grady,sz);
hessxx = reshape(hessxx,sz);
hessxy = reshape(hessxy,sz);
hessyy = reshape(hessyy,sz);
gradlapx = reshape(gradlapx,sz);
gradlapy = reshape(gradlapy,sz);

grad = cat(3,gradx,grady);
hess = cat(3,hessxx,hessxy,hessyy);
gradlap = cat(3,gradlapx,gradlapy);

end