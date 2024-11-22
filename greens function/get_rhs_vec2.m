function [rhs_vec] = get_rhs_vec2(coefs,kerns)

a0 = coefs{1};
abar = coefs{2}(:);
b0 = coefs{3};
bbar = coefs{4}(:);
g0 = coefs{5};
gbar = coefs{6}(:);
ax = coefs{7}(:);
ay = coefs{8}(:);
axx = coefs{9}(:);
axy = coefs{10}(:);
ayy = coefs{11}(:);
nu = coefs{end};

val = kerns{1};
hess = kerns{2};
gradlap = kerns{3};
phi = kerns{4};

hessxx = hess(:,:,1);
hessxy = hess(:,:,2);
hessyy = hess(:,:,3);

gradlapx = gradlap(:,:,1);
gradlapy = gradlap(:,:,2);

rhs = -(2*ax.*gradlapx + 2*ay.*gradlapy + (axx+ayy).*(hessxx+hessyy)+...
    (1-nu)*(2*axy.*hessxy - axx.*hessyy - ayy.*hessxx) + (abar/a0*b0-bbar).*val ) + ...
    - (gbar - g0*abar/a0).*phi;
rhs_vec = rhs(:);


end