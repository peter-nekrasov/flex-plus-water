function [val,grad,hess,third] = green(x,y,beta,gamma,surfeval)
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

[rts2, ejs] = find_roots(beta,gamma);

if surfeval
    ejs = ejs./rts2;
end

val = 0;

src = [0; 0];
targ = [r; 0*r];

for i = 1:5
    rhoj = rts2(i);
    ej = ejs(i);
    [ck0,~] = struveK(-rhoj, r);
    if rhoj == 0
        ck0 = 0;
    elseif angle(rhoj) == 0
       [ck0,~] = struveK(rhoj, r);
       [h0, ~] = helmdiffgreen(rhoj,src,targ);
       h0(r == 0) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));
       ck0 = -ck0 + 8*h0.' ; 
    end
    val = val + pi/2*ej*rhoj^2*ck0;
end

val = reshape(val,sz);

end