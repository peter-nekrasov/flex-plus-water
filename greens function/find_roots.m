function [rts,ejs] = find_roots(beta,gamma)
%%%%%
%
% Finds roots of polynomial z^5 - beta*z + gamma = 0
% 
% Returns:
% - rts: vector of roots of polynomial
% - ejs: coefficients in partial fraction expansion
%
%%%%%

    d1 = - beta;
    d0 = gamma;
    
    fp = @(z) 5*z.^4 + d1;
    
    C = bring_companion(d1,d0);
    
    rts1 = eig(C);
    rts = bring_refine_rts(rts1,d1,d0);
    
    ejs = 1./fp(rts); % partial fraction coefficients for 1/f

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