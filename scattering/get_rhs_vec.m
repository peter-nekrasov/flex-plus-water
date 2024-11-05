function [rhs_vec, rhs] = get_rhs_vec(coefs,k1,k2,phiinc)

k = sqrt(k1^2+k2^2);

abar = coefs{2};
bbar = coefs{4};
gbar = coefs{6};
ax = coefs{7};
ay = coefs{8};
axx = coefs{9};
axy = coefs{10};
ayy = coefs{11};
nu = coefs{end};

rhs = (-abar*k^5 + 2i*k^3*k1*ax + 2i*k^3*k2*ay + k^3*(axx+ayy)+ (1-nu)*k*(2*axy*k1*k2-axx*k2^2 - ayy*k1^2) + k*bbar - gbar).*phiinc;
rhs_vec = rhs(:);


end