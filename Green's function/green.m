function [val] = green(r,beta,gamma,surfeval)
%
% computes the green's function for the integro-differential equation
% whose solutions occur at the roots of the polynomial:
%             z^5 - beta*z + gamma = 0
%
% input:
%
% r - array of distances
% beta - coefficient beta in the equation
% gamma - coefficient gamma in the equation
%
% optional input:
%
% surfeval - boolean, default: false. If true, gives the kernel of the
%   Laplace single layer operator composed with the green's function needed
%   to evaluate the velocity potential phi on surface
%

if nargin < 4
    surfeval = false;
end

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
       ck0 = -ck0 + 8*h0.' ; 
    end
    val = val + pi/2*ej*rhoj^2*ck0;
end

val = reshape(val,sz);

end