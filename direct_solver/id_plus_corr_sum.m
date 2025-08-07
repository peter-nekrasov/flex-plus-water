function out = id_plus_corr_sum(V,corr,dinds,h)

a0 = V{1};
abar = V{2}(dinds);
b0 = V{3};
bbar = V{4}(dinds);
g0 = V{5};
alphax = V{7}(dinds);
alphay = V{8}(dinds);
alphaxx = V{9}(dinds);
alphaxy = V{10}(dinds);
alphayy = V{11}(dinds);
nu = V{12};

alphalap = alphaxx + alphayy;

Gs_hess_corr = corr{2};
Gs_gradlap_corr = corr{3};

Gs_corr = corr{1}(dinds,dinds)*h^2;
Gphi_corr = corr{4}(dinds,dinds)*h^2;

Gs_xx_corr = Gs_hess_corr{1}(dinds,dinds)*h^2;
Gs_xy_corr = Gs_hess_corr{2}(dinds,dinds)*h^2;
Gs_yy_corr = Gs_hess_corr{3}(dinds,dinds)*h^2;

Gslapx_corr = Gs_gradlap_corr{1}(dinds,dinds)*h^2;
Gslapy_corr = Gs_gradlap_corr{2}(dinds,dinds)*h^2;

Gslap_corr = Gs_xx_corr + Gs_yy_corr;

out = (a0 + abar)./a0.*speye(size(Gs_corr)) + alphax.*Gslapx_corr + alphay.*Gslapy_corr...
+ 0.5.*alphalap.*Gslap_corr + ...
+ 0.5*(1-nu)*(2*alphaxy.*Gs_xy_corr-alphayy.*Gs_xx_corr-alphaxx.*Gs_yy_corr) ...
- 0.5*(a0*bbar-abar*b0)./a0.*Gs_corr - 0.5*abar.*g0./a0.*Gphi_corr ;

end