%%%
%
% Testing my understanding by using zeta corrected trapezoid rule to
%       integrate 1/|r| against some compactly supported density
%
%%%

addpath('zetafunc/')

kern = @(x,y) 1./sqrt(x.^2 + y.^2); % kernel function
dens = @(x,y) exp(-(x.^2 + y.^2)/(2*5^2)); % test density 

trueint = integral2(@(x,y) kern(x,y).*dens(x,y),-50,50,-50,50,'AbsTol',0,'RelTol',1e-12);

% Trapezoid rule 

hs = 1./(2.^[1,2,3,4,5,6,7]);
errs0 = hs*0;
errs1 = hs*0;
z0 = epstein_zeta(1,1,0,1);

for ii = 1:numel(hs)
    h = hs(ii);
    [X,Y] = meshgrid(-50:h:50);
    kernmat = kern(X,Y).*dens(X,Y);
    kernmat(isnan(kernmat)) = 0;
    kernmat(isinf(kernmat)) = 0;
    errs0(ii) = abs(sum(kernmat(:)*h^2) - trueint);

    ind = intersect(find(X == 0),find(Y==0));
    [zi,zj] = ind2sub(size(X),ind);

    kernmat(zi,zj) = kernmat(zi,zj) - z0*dens(0,0)/h;
    errs1(ii) = abs(sum(kernmat(:)*h^2) - trueint);
end

figure(1)

loglog(hs,errs0,'-o')
hold on

loglog(hs,errs1,'-o')
hold on

loglog(hs,hs)
hold on

loglog(hs,0.01*hs.^3)

legend('without zeta correction','with zeta correction','h','h^3')

%% Now let's try to integrate log(|r|)

addpath('zetafunc/')

kern = @(x,y) log(sqrt(x.^2 + y.^2)); % kernel function
dens = @(x,y) exp(-(x.^2 + y.^2)/(2*5^2)); % test density 

trueint = integral2(@(x,y) kern(x,y).*dens(x,y),-40,40,-40,40,'AbsTol',0,'RelTol',1e-14);

% Trapezoid rule 

hs = 1./(2.^(1:8));
errs0 = hs*0;
errs1 = hs*0;
errs2 = hs*0;
z0 = imag(epstein_zeta(0+1i*10^-8,1,0,1))*1e8 ;

for ii = 1:numel(hs)

    h = hs(ii);
    z1 = log(h) ;

    [X,Y] = meshgrid(-40:h:40);
    kernmat = kern(X,Y).*dens(X,Y);

    ind = find((X == 0) & (Y==0));
    [idx,idy] = ind2sub(size(X),ind);

    kernmat(ind) = 0;
    errs0(ii) = abs(sum(kernmat(:)*h^2) - trueint);

    kernmat(ind) = kernmat(ind) + (z0+z1)*dens(0,0);
    errs1(ii) = abs(sum(kernmat(:)*h^2) - trueint);

    kernmat(ind)
end

figure(2)

loglog(hs,errs0,'-o')
hold on

loglog(hs,errs1,'-o')
hold on

loglog(hs,hs.^2)
hold on

loglog(hs,hs.^4)
hold on


legend('without zeta correction','with zeta correction','h^2','h^4')


rmpath('zetafunc/')
