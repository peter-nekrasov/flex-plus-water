function out = greenval(x,y,beta,gamma,opt)
%
% computes the green's function centered at (x,y) = 0 for the 
% integro-differential equation determined by the roots of the polynomial:
%             z^5 - beta*z + gamma = 0
%
% returns just the value of the Green's function (no derivatives)
% - out is the value of the Green's function centered at zero and
%   evaluated at (x,y)
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
% opt - bool, default: false. 
%         Possible options are:
%         opt = false => Green's function (and derivatives)
%         opt = true => kernel used to evaluate phi on surface
%

if nargin < 5
    opt = false;
end

r = sqrt(x.^2 + y.^2);
sz = size(r);
r = abs(r);
r = r(:).';
x = x(:).';
y = y(:).';

[rts2, ejs] = find_roots(beta,gamma);

if opt
    ejs = ejs./rts2;
end

val = 0;

src = [0; 0];
targ = [x; y];

for i = 1:5
    
    rhoj = rts2(i);
    ej = ejs(i);

    if angle(rhoj) == 0

       [sk0,~] = struveKdiffgreen(rhoj,src,targ);
       [h0,~] = helmdiffgreen(rhoj,src,targ);

       h0(r == 0) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));

       h0 = -4i*h0;

       val = val + ej*rhoj^2*(-sk0 + 2i*h0);

    elseif rhoj ~= 0

       [sk0,~] = struveKdiffgreen(-rhoj,src,targ);

       val = val + ej*rhoj^2*sk0;

    end

end

val = pi/2*val;

val = reshape(val,sz);

out = val;

end