function out = kern_sum(src,targ)

a0 = targ.a0;
b0 = targ.b0;
g0 = targ.g0;
abar = targ.abar;
alphax = targ.alphax;
alphay = targ.alphay;
alphaxx = targ.alphaxx;
alphaxy = targ.alphaxy;
alphayy = targ.alphayy;
alphalap = alphaxx+alphayy;
bbar = targ.bbar;
nu = targ.nu;

[rts,ejs] = find_roots(b0 / a0, g0 / a0);
ejs = ejs/a0;

kerns = green(src.r,targ.r,rts,ejs);

Gs = kerns{1};
Gs_hess = kerns{2};
Gs_gradlap = kerns{3};
Gphi = kerns{4};

Gsxx = Gs_hess(:,:,1);
Gsxy = Gs_hess(:,:,2);
Gsyy = Gs_hess(:,:,3);

Gslap = Gsxx+Gsyy;

Gslapx = Gs_gradlap(:,:,1);
Gslapy = Gs_gradlap(:,:,2);

out = alphax.*Gslapx + alphay.*Gslapy ...
+ 0.5.*alphalap.*Gslap + ...
+ 0.5*(1-nu)*(2*alphaxy.*Gsxy-alphayy.*Gsxx-alphaxx.*Gsyy) ...
- 0.5*(a0*bbar-abar*b0)./a0.*Gs - 0.5*abar.*g0./a0.*Gphi ;

end