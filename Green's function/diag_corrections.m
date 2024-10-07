h = 0.1;
xs = -10:h:10;
[~,n] = size(xs);
[X,Y] = meshgrid(xs);
beta = 5;
gamma = -1;

% Correction to Greens function
R = sqrt(X.^2 + Y.^2);
gval = green(R,beta,gamma);
gval(isnan(gval)) = 0;
testf = exp(-1/2*R.^2/(5*h^2));
integrand = @(x,y) green(sqrt(x.^2 + y.^2),beta,gamma).*exp(-1/2*(x.^2+y.^2)/(5*h^2));
beta_lg = integral2(integrand,min(xs),max(xs),min(xs),max(xs)) - sum(gval(:).*testf(:))*h^2;
gval(ceil(n/2),ceil(n/2)) = beta_lg/h^2;