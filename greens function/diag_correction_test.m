hs = [5,2,1,0.5,0.25,0.125,0.0625];
errs_no_correction = hs*0;
errs_w_correction = hs*0;
xmin = 0;
xmax = 50;
beta = 5;
gamma = -1;
integrand = @(x,y) green(sqrt(x.^2 + y.^2),beta,gamma,false).*exp(-1/2*(x.^2+y.^2));
true_int = integral2(integrand,xmin,xmax,xmin,xmax);

for ii = 1:numel(hs)
    h = hs(ii);
    xs = xmin:h:xmax;
    [~,n] = size(xs);
    [X,Y] = meshgrid(xs);

    % Correction to Greens function
    R = sqrt(X.^2 + Y.^2);
    gval = green(R,beta,gamma,false);
    testf = exp(-1/2*(R.^2));
    int_no_correction = sum(testf.*gval,'all').*h^2;
    errs_no_correction(ii) = abs(true_int - int_no_correction);

    gmat_corrected = center_correction(@(r) green(),h);
    %errs_w_correction(ii) = abs(true_int - int_w_correction);
end

loglog(hs, errs_no_correction)
hold on
loglog(hs, hs.^4)