%%%
%
% Testing my understanding by using zeta corrected trapezoid rule to
%       integrate 1/|r| against some compactly supported density
%
%%%

addpath('zetafunc/')

kern = @(x,y) 1./sqrt(x.^2 + y.^2); % kernel function
dens = @(x,y) exp(-25*((x).^2 + (y).^2).^4 - 5*x.^2); % test density 

trueint = integral2(@(x,y) dens(x,y).*kern(x,y),-1,1,-1,1,"AbsTol",0,"RelTol",1e-16)

%trueint = 1.8725427747132428474951878691780619935404005243656784679539306;

% Trapezoid rule 

hs = 1./(2.^(1:11));
errs0 = hs*0;
errs1 = hs*0;
errs2 = hs*0;
errs3 = hs*0;
[z0] = epstein_zeta_int(1,1,0,1);
[~,z1] = epstein_zeta_int(-1,1,0,1,1,0,0);
[~,~,z2] = epstein_zeta_int(-3,1,0,1,1,0,0);
[~,~,z3] = epstein_zeta_int(-3,1,0,1,0,1,0);

for ii = 1:numel(hs)
    h = hs(ii);
    [X,Y] = meshgrid(-1:h:1);
    kernmat = kern(X,Y)*h^2;
    kernmat(isnan(kernmat)) = 0;
    kernmat(isinf(kernmat)) = 0;
    errs0(ii) = abs(sum(kernmat.*dens(X,Y),'all') - trueint)/abs(trueint);

    ind = intersect(find(X == 0),find(Y==0));
    [zi,zj] = ind2sub(size(X),ind);

    kernmat(zi,zj) = kernmat(zi,zj) - z0*dens(0,0)*h;
    errs1(ii) = abs(sum(kernmat.*dens(X,Y),'all') - trueint)/abs(trueint);

    kernmat = kern(X,Y)*h^2;
    kernmat(isnan(kernmat)) = 0;
    kernmat(isinf(kernmat)) = 0;


    A = [1 1 1 1 1;
        0 1 -1 0 0;
        0 0 0 1 -1;
        0 1 1 0 0;
        0 0 0 1 1];
    b = [-h*z0; 0; 0; -2*h*z1; -2*h*z1];
    tau = A \ b;

    kernmat(zi,zj) = kernmat(zi,zj) + dens(0,0)*tau(1);
    kernmat(zi+1,zj) = kernmat(zi+1,zj) + tau(2);
    kernmat(zi-1,zj) = kernmat(zi-1,zj) + tau(3);
    kernmat(zi,zj+1) = kernmat(zi,zj+1) + tau(4);
    kernmat(zi,zj-1) = kernmat(zi,zj-1) + tau(5);
    errs2(ii) = abs(sum(kernmat.*dens(X,Y),'all') - trueint)/trueint;

    kernmat = kern(X,Y)*h^2;
    kernmat(isnan(kernmat)) = 0;
    kernmat(isinf(kernmat)) = 0;

    A = [ones(1,13); ...
        0 1 -1 0 0 2 -2 0 0 1 -1 1 -1; ...
        0 0 0 1 -1 0 0 2 -2 1 1 -1 -1; ...
        0 1 1 0 0 4 4 0 0 1 1 1 1; ...
        0 0 0 1 1 0 0 4 4 1 1 1 1; ...
        zeros(1,9) 1 -1 -1 1; ...
        zeros(1,9) 1 1 -1 -1; ...
        zeros(1,9) 1 -1 1 -1; ...
        0 1 -1 0 0 8 -8 0 0 1 -1 1 -1; ...
        0 0 0 1 -1 0 0 8 -8 1 1 -1 -1; ...
        0 1 1 0 0 16 16 0 0 1 1 1 1; ...
        zeros(1,9) ones(1,4); ...
        0 0 0 1 1 0 0 16 16 1 1 1 1];
    b = [z0; 0; 0; 2*z1; 2*z1; 0; 0; 0; 0; 0; 4/3*z2; 1/3*z3; 4/3*z2];
    b = -h*b;
    tau = A \ b;

    kernmat(zi,zj) = kernmat(zi,zj) + tau(1);
    kernmat(zi+1,zj) = kernmat(zi+1,zj) + tau(2);
    kernmat(zi-1,zj) = kernmat(zi-1,zj) + tau(3);
    kernmat(zi,zj+1) = kernmat(zi,zj+1) + tau(4);
    kernmat(zi,zj-1) = kernmat(zi,zj-1) + tau(5);
    kernmat(zi+2,zj) = kernmat(zi+2,zj) + tau(6);
    kernmat(zi-2,zj) = kernmat(zi-2,zj) + tau(7);
    kernmat(zi,zj+2) = kernmat(zi,zj+2) + tau(8);
    kernmat(zi,zj-2) = kernmat(zi,zj-2) + tau(9);
    kernmat(zi+1,zj+1) = kernmat(zi+1,zj+1) + tau(10);
    kernmat(zi-1,zj+1) = kernmat(zi-1,zj+1) + tau(11);
    kernmat(zi+1,zj-1) = kernmat(zi+1,zj-1) + tau(12);    
    kernmat(zi-1,zj-1) = kernmat(zi-1,zj-1) + tau(13);

    errs3(ii) = abs(sum(kernmat.*dens(X,Y),'all') - trueint)/trueint;

end

figure(1)

loglog(hs,errs0,'-o')
hold on

loglog(hs,errs1,'-o')
hold on

loglog(hs,errs2,'-o')
hold on

loglog(hs,errs3,'-o')
hold on

loglog(hs,0.5*hs,'--')
hold on

loglog(hs,0.1*hs.^3,'--')
hold on

loglog(hs,0.3*hs.^5,'--')
hold on

loglog(hs,0.1*hs.^7,'--')
hold on

xlim([min(hs) max(hs)])
ylim([0.1*min(errs3) 10*max(errs3)])

legend('w/o zeta correction','1 pt zeta correction', '5 pt zeta correction', '13 pt zeta correction', 'h','h^3','h^5','h^7','Location','best')
title('zeta correction for 1/r')


%% Now let's try to integrate log(|r|)

addpath('zetafunc/')

kern = @(x,y) log(sqrt(x.^2 + y.^2)); % kernel function
dens = @(x,y) exp(-25*((x).^2 + (y).^2).^4 - 5*x.^2); % test density 

trueint = integral2(@(x,y) dens(x,y).*kern(x,y),-1,1,-1,1,"AbsTol",0,"RelTol",1e-16)

% Trapezoid rule 

hs = 1./(2.^(0:8));
errs0 = hs*0;
errs1 = hs*0;
errs2 = hs*0;
errs3 = hs*0;
[z0] = imag(epstein_zeta(0+1i*10^-12,1,0,1))*1e12 ;
[~,z1] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
z1 = imag(z1)*1e12;
[~,~,z2] = epstein_zeta(-4+1i*10^-12,1,0,1,1,0,0) ;
z2 = imag(z2)*1e12;
[~,~,z3] = epstein_zeta(-4+1i*10^-12,1,0,1,0,1,0) ;
z3 = imag(z3)*1e12;

for ii = 1:numel(hs)

    h = hs(ii);

    [X,Y] = meshgrid(-40:h:40);
    kernmat = kern(X,Y)*h^2;

    ind = find((X == 0) & (Y==0));
    [zi,zj] = ind2sub(size(X),ind);

    kernmat(ind) = 0;
    errs0(ii) = abs(sum(kernmat.*dens(X,Y),'all') - trueint)/abs(trueint);

    kernmat(ind) = kernmat(ind) + (z0+log(h))*h^2;
    errs1(ii) = abs(sum(kernmat.*dens(X,Y),'all') - trueint)/abs(trueint);

    kernmat = kern(X,Y)*h^2;
    kernmat(ind) = 0;

    A = [1 1 1 1 1;
        0 1 -1 0 0;
        0 0 0 1 -1;
        0 1 1 0 0;
        0 0 0 1 1];
    b = [h^2*z0 + h^2*log(h); 0; 0; h^2*z1; h^2*z1];
    tau = A \ b;

    kernmat(zi,zj) = kernmat(zi,zj) + tau(1);
    kernmat(zi+1,zj) = kernmat(zi+1,zj) + tau(2);
    kernmat(zi-1,zj) = kernmat(zi-1,zj) + tau(3);
    kernmat(zi,zj+1) = kernmat(zi,zj+1) + tau(4);
    kernmat(zi,zj-1) = kernmat(zi,zj-1) + tau(5);
    errs2(ii) = abs(sum(kernmat.*dens(X,Y),'all') - trueint)/abs(trueint);

    kernmat = kern(X,Y)*h^2;
    kernmat(ind) = 0;

    A = [ones(1,13); ...
        0 1 -1 0 0 2 -2 0 0 1 -1 1 -1; ...
        0 0 0 1 -1 0 0 2 -2 1 1 -1 -1; ...
        0 1 1 0 0 4 4 0 0 1 1 1 1; ...
        0 0 0 1 1 0 0 4 4 1 1 1 1; ...
        zeros(1,9) 1 -1 -1 1; ...
        zeros(1,9) 1 1 -1 -1; ...
        zeros(1,9) 1 -1 1 -1; ...
        0 1 -1 0 0 8 -8 0 0 1 -1 1 -1; ...
        0 0 0 1 -1 0 0 8 -8 1 1 -1 -1; ...
        0 1 1 0 0 16 16 0 0 1 1 1 1; ...
        zeros(1,9) ones(1,4); ...
        0 0 0 1 1 0 0 16 16 1 1 1 1];
    b = [z0 + log(h); 0; 0; z1; z1; 0; 0; 0; 0; 0; 1/2*z2; 1/8*z3; 1/2*z2];
    b = h^2*b;
    tau = A \ b;

    kernmat(zi,zj) = kernmat(zi,zj) + tau(1);
    kernmat(zi+1,zj) = kernmat(zi+1,zj) + tau(2);
    kernmat(zi-1,zj) = kernmat(zi-1,zj) + tau(3);
    kernmat(zi,zj+1) = kernmat(zi,zj+1) + tau(4);
    kernmat(zi,zj-1) = kernmat(zi,zj-1) + tau(5);
    kernmat(zi+2,zj) = kernmat(zi+2,zj) + tau(6);
    kernmat(zi-2,zj) = kernmat(zi-2,zj) + tau(7);
    kernmat(zi,zj+2) = kernmat(zi,zj+2) + tau(8);
    kernmat(zi,zj-2) = kernmat(zi,zj-2) + tau(9);
    kernmat(zi+1,zj+1) = kernmat(zi+1,zj+1) + tau(10);
    kernmat(zi-1,zj+1) = kernmat(zi-1,zj+1) + tau(11);
    kernmat(zi+1,zj-1) = kernmat(zi+1,zj-1) + tau(12);    
    kernmat(zi-1,zj-1) = kernmat(zi-1,zj-1) + tau(13);

    errs3(ii) = abs(sum(kernmat.*dens(X,Y),'all') - trueint)/abs(trueint);

end

figure(2)

loglog(hs,errs0,'-o')
hold on

loglog(hs,errs1,'-o')
hold on

loglog(hs,errs2,'-o')
hold on

loglog(hs,errs3,'-o')
hold on

loglog(hs,hs.^2,'--')
hold on

loglog(hs,hs.^4,'--')
hold on

loglog(hs,0.5*hs.^6,'--')
hold on

loglog(hs,1.3*hs.^8,'--')
hold on

xlim([min(hs) max(hs)])
ylim([0.1*min(errs2) 10*max(errs2)])


legend('w/o zeta correction','1 pt correction','5 pt correction','13 pt correction','h^2','h^4','h^6','h^8','Location','southeast')
title('zeta correction for log(|r|)')


rmpath('zetafunc/')


%% Now let's try to integrate \vec{r} / r^2 !!!!!

addpath('zetafunc/')

kern = @(x,y) x./(x.^2 + y.^2); % kernel function
dens = @(x,y) (cos(x)+sin(x)).*exp(-30*((x).^2 + (y).^2)); % test density 

trueint = integral2(@(x,y) dens(x,y).*kern(x,y),-1,1,-1,1,"AbsTol",0,"RelTol",1e-16)

% Trapezoid rule 

hs = 1./(2.^(0:8));
errs0 = hs*0;
errs1 = hs*0;
errs2 = hs*0;
errs3 = hs*0;

[~,z1] = epstein_zeta(0+1i*10^-12,1,0,1,1,0,0) ;
z1 = imag(z1)*1e12;

for ii = 1:numel(hs)

    h = hs(ii);

    [X,Y] = meshgrid(-40:h:40);
    kernmat = kern(X,Y)*h^2;

    ind = find((X == 0) & (Y==0));
    [zi,zj] = ind2sub(size(X),ind);

    kernmat(ind) = 0;
    errs0(ii) = abs(sum(kernmat.*dens(X,Y),'all') - trueint)/abs(trueint);

    kernmat = kern(X,Y)*h^2; 
    kernmat(ind) = 0;

    A = [1 1 1 1 1;
        0 1 -1 0 0;
        0 0 0 1 -1;
        0 1 1 0 0;
        0 0 0 1 1];
    b = [0; 2*z1; 0; 0; 0];
    b = b*h;
    tau = A \ b;

    kernmat(zi,zj) = kernmat(zi,zj) + tau(1);
    kernmat(zi+1,zj) = kernmat(zi+1,zj) + tau(2);
    kernmat(zi-1,zj) = kernmat(zi-1,zj) + tau(3);
    kernmat(zi,zj+1) = kernmat(zi,zj+1) + tau(4);
    kernmat(zi,zj-1) = kernmat(zi,zj-1) + tau(5);
    errs1(ii) = abs(sum(kernmat.*dens(X,Y),'all') - trueint)/abs(trueint);

    kernmat = kern(X,Y)*h^2;
    kernmat(ind) = 0;

    % 
    % A = [ones(1,13); ...
    %     0 1 -1 0 0 2 -2 0 0 1 -1 1 -1; ...
    %     0 0 0 1 -1 0 0 2 -2 1 1 -1 -1; ...
    %     0 1 1 0 0 4 4 0 0 1 1 1 1; ...
    %     0 0 0 1 1 0 0 4 4 1 1 1 1; ...
    %     zeros(1,9) 1 -1 -1 1; ...
    %     zeros(1,9) 1 1 -1 -1; ...
    %     zeros(1,9) 1 -1 1 -1; ...
    %     0 1 -1 0 0 8 -8 0 0 1 -1 1 -1; ...
    %     0 0 0 1 -1 0 0 8 -8 1 1 -1 -1; ...
    %     0 1 1 0 0 16 16 0 0 1 1 1 1; ...
    %     zeros(1,9) ones(1,4); ...
    %     0 0 0 1 1 0 0 16 16 1 1 1 1];
    % b = [z0 + log(h); 0; 0; z1; z1; 0; 0; 0; 0; 0; 1/2*z2; 1/8*z3; 1/2*z2];
    % b = h^2*b;
    % tau = A \ b;
    % 
    % kernmat(zi,zj) = kernmat(zi,zj) + tau(1);
    % kernmat(zi+1,zj) = kernmat(zi+1,zj) + tau(2);
    % kernmat(zi-1,zj) = kernmat(zi-1,zj) + tau(3);
    % kernmat(zi,zj+1) = kernmat(zi,zj+1) + tau(4);
    % kernmat(zi,zj-1) = kernmat(zi,zj-1) + tau(5);
    % kernmat(zi+2,zj) = kernmat(zi+2,zj) + tau(6);
    % kernmat(zi-2,zj) = kernmat(zi-2,zj) + tau(7);
    % kernmat(zi,zj+2) = kernmat(zi,zj+2) + tau(8);
    % kernmat(zi,zj-2) = kernmat(zi,zj-2) + tau(9);
    % kernmat(zi+1,zj+1) = kernmat(zi+1,zj+1) + tau(10);
    % kernmat(zi-1,zj+1) = kernmat(zi-1,zj+1) + tau(11);
    % kernmat(zi+1,zj-1) = kernmat(zi+1,zj-1) + tau(12);    
    % kernmat(zi-1,zj-1) = kernmat(zi-1,zj-1) + tau(13);
    % 
    % errs3(ii) = abs(sum(kernmat.*dens(X,Y),'all') - trueint)/abs(trueint);

end

figure(2)

loglog(hs,errs0,'-o')
hold on

loglog(hs,errs1,'-o')
hold on
% 
% loglog(hs,errs2,'-o')
% hold on
% 
% loglog(hs,errs3,'-o')
% hold on
% 
loglog(hs,5*hs.^2,'--')
hold on
% 
% loglog(hs,hs.^4,'--')
% hold on
% 
% loglog(hs,0.5*hs.^6,'--')
% hold on
% 
% loglog(hs,1.3*hs.^8,'--')
% hold on
% 
% xlim([min(hs) max(hs)])
% ylim([0.1*min(errs2) 10*max(errs2)])


legend('w/o zeta correction','1 pt correction','5 pt correction','13 pt correction','h^2','h^4','h^6','h^8','Location','southeast')
title('zeta correction for x/r^2')


rmpath('zetafunc/')