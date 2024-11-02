%% Checking corrections for G 

addpath(genpath('..'))

gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);
targ = [2; 2];

dens = @(x,y) x.*exp(-(x.^2+y.^2)/(5));
greenfac = @(x,y) greenvalonly(targ,x,y,rts,ejs);
truev = -0.305656416766226 - 0.016380210577180i; % integral2(@(x,y) dens(x,y).*greenfac(x,y),-20,20,-20,20,"AbsTol",0,"RelTol",10E-16);

hs = [2 1 0.5 0.25 0.2 0.1 0.05];
errs0 = hs*0;
errs1 = hs*0;

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-12:h:12);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs));
    val = kern{1};
    d1 = dens(X,Y);
    dint = sum(val.'.*d1(:));
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(rts,ejs,h);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),inds,corrs);
    val = kern{1};
    dint = sum(val.*d1(:).');
    errs1(ii) = abs(dint - truev);

end

loglog(hs,errs0,'o-')
hold on

loglog(hs,errs1,'o-')
hold on

loglog(hs,0.0005*hs.^4,'--')
hold on

loglog(hs,0.0001*hs.^6,'--')
hold on 

legend('no correction', '5 pt correction', 'h^4', 'h^6') 

%% Checking corrections for G_{xx} 

addpath(genpath('..'))

gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);
targ = [2; 2];

dens = @(x,y) exp(-(x.^2+y.^2)/(10));
truev = -0.002078708423374 + 0.012091284743428i;
%integral2(@(x,y) dens(x,y).*hessxxonly(targ,x,y,rts,ejs),-15,15,-15,15,'AbsTol',0,'RelTol',10E-14)

hs = [2 1 0.5 0.2 0.1 0.05 0.025];
errs0 = hs*0;
errs1 = hs*0; 

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-15:h:15);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs));
    val = kern{2};
    val = val(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(rts,ejs,h);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),inds,corrs);
    val = kern{2};
    val = val(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev);

end


loglog(hs,errs0,'o-')
hold on

loglog(hs,errs1,'o-')
hold on

loglog(hs,0.1*hs.^2,'--')
hold on

loglog(hs,0.0005*hs.^6,'--')
hold on

legend('no correction', '5 pt correction', 'h^2', 'h^6','Location','northwest') 


%% Checking corrections for G_{xy} 

addpath(genpath('..'))

gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);
targ = [2; 2];

dens = @(x,y) sin(x+y).*exp(-(x.^2+y.^2)/(10));
truev = 0.275566788162816 + 0.239570256490271i; %integral2(@(x,y) dens(x,y).*hessxyonly(targ,x,y,rts,ejs),-15,15,-15,15,'AbsTol',0,'RelTol',10E-14)

hs = [2 1 0.5 0.2 0.1 0.05 0.025];
errs0 = hs*0;
errs1 = hs*0; 

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-15:h:15);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs));
    val = kern{2};
    val = val(:,:,2);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(rts,ejs,h);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),inds,corrs);
    val = kern{2};
    val = val(:,:,2);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev);

end

loglog(hs,errs0,'o-')
hold on

loglog(hs,errs1,'o-')
hold on

loglog(hs,0.1*hs.^4,'--')
hold on

loglog(hs,0.0005*hs.^6,'--')
hold on

legend('no correction', '13 pt correction', 'h^2', 'h^6','Location','northwest') 

%% Checking corrections for G_{yy} 

addpath(genpath('..'))

gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);
targ = [2; 2];

dens = @(x,y) exp(-(x.^2+y.^2)/(10));
%truev =   0.057816396686963 - 0.015079981170599i;
truev = -0.002078708423374 + 0.012091284743428i;
%integral2(@(x,y) dens(x,y).*hessxxonly(targ,x,y,rts,ejs),-15,15,-15,15,'AbsTol',0,'RelTol',10E-14)

hs = [2 1 0.5 0.2 0.1 0.05 0.025];
errs0 = hs*0;
errs1 = hs*0; 

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-15:h:15);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs));
    val = kern{2};
    val = val(:,:,3);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(rts,ejs,h);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),inds,corrs);
    val = kern{2};
    val = val(:,:,3);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev);

end


loglog(hs,errs0,'o-')
hold on

loglog(hs,errs1,'o-')
hold on

loglog(hs,0.1*hs.^2,'--')
hold on

loglog(hs,0.0005*hs.^6,'--')
hold on

legend('no correction', '5 pt correction', 'h^2', 'h^6','Location','northwest') 

%%  Checking corrections for phi 

addpath(genpath('..'))

gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);
targ = [2; 2];

dens = @(x,y) exp(-(x.^2+y.^2)/(5));
greenfac = @(x,y) phivalonly(targ,x,y,rts,ejs);
%truev = -0.062850332948632 + 0.069122770940408i; %  integral2(@(x,y) dens(x,y).*greenfac(x,y),-50,50,-50,50,"AbsTol",0,"RelTol",10E-14);
truev = -0.236168737742442 - 0.050305328873419i; % integral2(@(x,y) dens(x,y).*greenfac(x,y),-20,20,-20,20,"AbsTol",0,"RelTol",10E-16);

hs = [2 1 0.5 0.25 0.2 0.1];
errs0 = hs*0;
errs1 = hs*0;

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-12:h:12);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs));
    val = kern{4};
    d1 = dens(X,Y);
    dint = sum(val.'.*d1(:));
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(rts,ejs,h);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),inds,corrs);
    val = kern{4};
    dint = sum(val.'.*d1(:));
    errs1(ii) = abs(dint - truev);

end

loglog(hs,errs0,'o-')
hold on

loglog(hs,errs1,'o-')
hold on

loglog(hs,0.0005*hs.^5,'--')
hold on

loglog(hs,0.0001*hs.^6,'--')
hold on 

legend('no correction', '5 pt correction', 'h^5', 'h^6') 


function val = greenvalonly(targ,x,y,rts,ejs)

    src = [x(:).'; y(:).'];
    gf = green(src,targ,rts,ejs);
    val = reshape(gf{1},size(x));
    

end

function val = hessxyonly(targ,x,y,rts,ejs)

    src = [x(:).'; y(:).'];
    gf = green(src,targ,rts,ejs);
    gf = gf{2};
    val = gf(:,:,2);
    val = reshape(val,size(x));

end

function val = hessxxonly(targ,x,y,rts,ejs)

    src = [x(:).'; y(:).'];
    gf = green(src,targ,rts,ejs);
    gf = gf{2};
    val = gf(:,:,1);
    val = reshape(val,size(x));

end

function val = phivalonly(targ,x,y,rts,ejs)

    src = [x(:).'; y(:).'];
    gf = green(src,targ,rts,ejs);
    val = reshape(gf{4},size(x));
    

end
