function [xmat] = kernbyindex(i,j,rs,h,corrs,rts,ejs,V)

src  = rs(:,j);
targ = rs(:,i);
vals = green(src,targ,rts,ejs);

v1 = vals{1} + corrs{1}(i,j);
v2a= squeeze(vals{2}(:,:,1)) + corrs{2}{1}(i,j);
v2b= squeeze(vals{2}(:,:,2)) + corrs{2}{2}(i,j);
v2c= squeeze(vals{2}(:,:,3)) + corrs{2}{3}(i,j);
v3a= squeeze(vals{3}(:,:,1)) + corrs{3}{1}(i,j);
v3b= squeeze(vals{3}(:,:,2)) + corrs{3}{2}(i,j);

v4 = vals{4} + corrs{4}(i,j);

v2lap = v2a + v2c;

%%%%%%%

a0 = V{1};
abar = V{2}(i);
b0 = V{3};
bbar = V{4}(i);
g0 = V{5};
gbar = V{6};
alphax = V{7}(i);
alphay = V{8}(i);
alphaxx = V{9}(i);
alphaxy = V{10}(i);
alphayy = V{11}(i);
nu = V{end};
alphalap = alphaxx + alphayy;


xmat = diag(alphax)*v3a + diag(alphay)*v3b + ...
    0.5*diag(alphalap)*v2lap + ...
    0.5*(1-nu)*(2*diag(alphaxy)*v2b-diag(alphayy)*v2a-diag(alphaxx)*v2c)-...
    0.5*diag((a0*bbar-abar*b0)./a0)*v1 - diag(0.5*abar.*g0./a0)*v4;
xmat = xmat*h^2;

[ind1,ind2] = find(i(:)==j(:).');
ind0 = find(i(:)==j(:).');

xmat(ind0) = xmat(ind0) + (a0 + abar(ind1).')./a0;

end