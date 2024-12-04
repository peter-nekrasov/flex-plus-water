addpath(genpath('../quadrature'))

h = 0.5;
[X, Y] = meshgrid(-10:h:10);
a0 = 2;

[inds,corrs] = get_correct(h,a0);

corrs = get_sparse_corr(size(X),inds,corrs);

pts = [X(:).'; Y(:).'];

b0 = 3;
g0 = -2;

[rts,ejs] = find_roots(b0 / a0, g0 / a0);
ejs = ejs/a0;

mats = green(pts,pts,rts,ejs);


return




