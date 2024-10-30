%% Checking corrections for G 

addpath(genpath('..'))

gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);
targ = [2; 2];

dens = @(x,y) x.*exp(-(x.^2+y.^2)/(5));
greenfac = @(x,y) greenvalonly(targ,x,y,rts,ejs);
%truev = -0.062850332948632 + 0.069122770940408i; %  integral2(@(x,y) dens(x,y).*greenfac(x,y),-50,50,-50,50,"AbsTol",0,"RelTol",10E-14);
truev = -0.152828208383113 - 0.008190105288590i; % integral2(@(x,y) dens(x,y).*greenfac(x,y),-20,20,-20,20,"AbsTol",0,"RelTol",10E-16);

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

dens = @(x,y) x.*exp(-(x.^2+y.^2)/(5));
greenfac = @(x,y) hessxxonly(targ,x,y,rts,ejs);
%truev = -0.062850332948632 + 0.069122770940408i; %  integral2(@(x,y) dens(x,y).*greenfac(x,y),-50,50,-50,50,"AbsTol",0,"RelTol",10E-14);
truev = integral2(@(x,y) dens(x,y).*greenfac(x,y),-12,12,-12,12,"AbsTol",0,"RelTol",10E-11);

hs = [2 1 0.5 0.25 0.2 0.1];
errs0 = hs*0;
errs1 = hs*0;

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-12:h:12);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs));
    val = kern{2};
    val = val(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val.'.*d1(:));
    errs0(ii) = abs(dint - truev);

    % [inds, corrs] = get_correct(rts,ejs,h);
    % kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),inds,corrs);
    % val = kern{2};
    % val = val(:,:,2)
    % dint = sum(val.'.*d1(:));
    % errs1(ii) = abs(dint - truev);

end

loglog(hs,errs0,'o-')
hold on

% loglog(hs,errs1,'o-')
% hold on

loglog(hs,0.0005*hs.^4,'--')
hold on

% loglog(hs,0.0001*hs.^6,'--')
% hold on 

legend('no correction', '5 pt correction', 'h^4', 'h^6','Location','northwest') 


%%  Checking corrections for phi 

addpath(genpath('..'))

gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);
targ = [2; 2];

dens = @(x,y) x.*exp(-(x.^2+y.^2)/(5));
greenfac = @(x,y) phivalonly(targ,x,y,rts,ejs);
%truev = -0.062850332948632 + 0.069122770940408i; %  integral2(@(x,y) dens(x,y).*greenfac(x,y),-50,50,-50,50,"AbsTol",0,"RelTol",10E-14);
truev = -0.186834150038740 - 0.005897287268752i; % integral2(@(x,y) dens(x,y).*greenfac(x,y),-20,20,-20,20,"AbsTol",0,"RelTol",10E-16);

hs = [2 1 0.5 0.25 0.2 0.1 0.05];
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

loglog(hs,0.0005*hs.^4,'--')
hold on

loglog(hs,0.0001*hs.^6,'--')
hold on 

legend('no correction', '5 pt correction', 'h^4', 'h^6') 


function val = greenvalonly(targ,x,y,rts,ejs)

    src = [x(:).'; y(:).'];
    gf = green(src,targ,rts,ejs);
    val = reshape(gf{1},size(x));
    

end

function val = hessxxonly(targ,x,y,rts,ejs)

    src = [x(:).'; y(:).'];
    gf = green(src,targ,rts,ejs);
    gf = gf{2};
    gf = gf(:,:,1);
    val = reshape(gf,size(x));
    

end

function val = phivalonly(targ,x,y,rts,ejs)

    src = [x(:).'; y(:).'];
    gf = green(src,targ,rts,ejs);
    val = reshape(gf{4},size(x));
    

end
