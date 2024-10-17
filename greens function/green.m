function [val,grad] = green(x,y,beta,gamma,surfeval)
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
% - third is the gradient of the Laplacian, namely 
%   third(:,:,1) has G_{xxx} + G_{xyy}, third(:,:,2) has G_{yxx} + G_{yyy}
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
% surfeval - boolean, default: false. If true, gives the kernel of the
%   Laplace single layer operator composed with the green's function needed
%   to evaluate the velocity potential phi on surface
%

if nargin < 5
    surfeval = false;
end

r = sqrt(x.^2 + y.^2);
sz = size(r);
r = abs(r);
r = r(:).';
x = x(:).';
y = y(:).';

[rts2, ejs] = find_roots(beta,gamma);

if surfeval
    ejs = ejs./rts2;
end

val = 0;
gradx = 0;
grady = 0;

src = [0; 0];
targ = [x; y];

for i = 1:5
    
    rhoj = rts2(i);
    ej = ejs(i);

    if angle(rhoj) == 0

       [ck0,ck1] = struveK(rhoj, r);
       [h0,gradh0] = helmdiffgreen(rhoj,src,targ);
       h0(r == 0) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));
       h0 = -4i*h0.';
       gradh0 = -4i*gradh0;
       h0x = gradh0(:,:,1).';
       h0y = gradh0(:,:,2).';

       val = val + ej*rhoj^2*(-ck0 + 2i*h0);

       gradx = gradx + ej*rhoj^3*x./r.*ck1 + ej*rhoj^2*(2i*h0x);
       grady = grady + ej*rhoj^3*y./r.*ck1 + ej*rhoj^2*(2i*h0y);

    elseif rhoj ~= 0

       [ck0,ck1] = struveK(-rhoj, r);
       
       val = val + ej*rhoj^2*ck0;

       gradx = gradx + ej*rhoj^3*x./r.*ck1;
       grady = grady + ej*rhoj^3*y./r.*ck1;

    end

end

val = pi/2*val;
gradx = pi/2*gradx;
grady = pi/2*grady;

val = reshape(val,sz);
gradx = reshape(gradx,sz);
grady = reshape(grady,sz);
grad = cat(3,gradx,grady);

end