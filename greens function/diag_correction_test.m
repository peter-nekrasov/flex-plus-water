hs = [4 2 1 0.5 0.25 0.125 0.125/2 0.125/4 0.125/8 0.125/16];
errs = hs*0;
xmin = -50;
xmax = 50;
beta = 5;
gamma = -1;
integrand = @(x,y) green(sqrt(x.^2 + y.^2),beta,gamma,false).*exp(-1/pi*(x.^2+y.^2));
true_int = integral2(integrand,xmin,xmax,xmin,xmax,'RelTol',1e-13,'AbsTol',0);

for ii = 1:numel(hs)
    h = hs(ii);
    xs = xmin:h:xmax;
    [~,n] = size(xs);
    [X,Y] = meshgrid(xs);

    % Correction to Greens function
    R = sqrt(X.^2 + Y.^2);
    gval = green(R,beta,gamma,false);
    testf = exp(-1/pi*(R.^2));
    intmat = kernmat(@(y1,y2) green(sqrt(y1.^2 + y2.^2),beta,gamma),h,xmax,1);
    int = sum(testf.*intmat,'all');
    errs(ii) = abs(true_int - int);
    disp(abs(true_int - int))

    %intmat = kernmat(@(x,y) green(sqrt(x.^2+y.^2),beta,gamma),h,[xmin xmax],[xmin xmax],5);
    %errs_w_correction(ii) = abs(true_int - int_corrected);
end

%% Plotting error

loglog(hs, errs,'o-')
hold on
loglog(hs, hs.^4)
legend('1 pt stencil','h^4')
ylabel('absolute error')
xlabel('h')

%% Plotting 

% loglog(hs, errs,'o-', hs, pts1_errs, 'o-', hs, pts5_errs, 'o-', hs, pts13_errs, 'o-')
% hold on
% loglog(hs,hs.^4,'--',hs,0.00005*hs.^10,'--')
% ylim([min(pts13_errs),max(pts13_errs)])
% legend('no correction', '1 pt stencil','5 pt stencil', '13 pt stencil','h^4', 'h^{10}' ,'Location','northwest')
% ylabel('relative error')
% xlabel('h')
