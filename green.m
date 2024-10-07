function [val] = green(r,beta,gamma)
%
% computing the green's function for  our pde
% whose solutions occur at the roots of the polynomial:
%             z^5 - beta*z + gamma = 0
%

sz = size(r);
r = abs(r);
r = r(:).';

d1 = - beta;
d0 = gamma;

f = @(z) z.^5 + d1*z + d0;
fp = @(z) 5*z.^4 + d1;

C = bring_companion(d1,d0);

rts1 = eig(C);
rts2 = bring_refine_rts(rts1,d1,d0);

ejs = 1./fp(rts2); % partial fraction coefficients for 1/f

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

function C = bring_companion(d1,d0)
%
% returns the companion matrix of z^5 + d1 z + d0
%

C = diag(ones(4,1),-1);
C(:,end) = [-d0; -d1; 0; 0; 0];

end

function rts2 = bring_refine_rts(rts1,d1,d0,nstep)
%
% applies nstep steps of Newton to each root
%

if nargin < 4
    nstep = 3;
end

f = @(z) z.^5 + d1*z + d0;
fp = @(z) 5*z.^4 + d1;

rts2 = rts1;

for j = 1:nstep
    rts2 = rts2 - f(rts2)./fp(rts2);
end

end