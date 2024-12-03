h = 0.01;
[X, Y] = meshgrid(-10:h:10);
a0 = 2;

[inds,corrs] = get_correct(h,a0);

corrs = get_sparse_corr(size(X),inds,corrs);

return




