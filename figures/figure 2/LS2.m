%% Checking corrections for \partial_x (\Delta G) 

% addpath(genpath('..'))

gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);
targ = [2; 2];

dens = @(x,y) x.*exp(-(x.^2+y.^2)/2).*cos(y/10+1);
%truev =   0.057816396686963 - 0.015079981170599i;
% truev = 1.027910069045180 + 4.450033610594533i; % 
truev = integral2(@(x,y) dens(x,y).*gradlapxonly(targ,x,y,rts,ejs),-10,10,-10,10,'AbsTol',10E-18,'RelTol',10E-18);

hs = [2 1 0.5 0.2 0.1 0.05 0.025 0.025/2 0.025/4 0.025/8];
errs0 = hs*0;
errs1 = hs*0; 

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-10:h:10);
    src = [X(:).'; Y(:).'];
    % kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
    % val = kern{3};
    % val = val(:,:,1);
    % d1 = dens(X,Y);
    % dint = sum(val(:).*d1(:),'all');
    % errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(h,1);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
    val = kern{3};
    val = val(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end
K1 = errs1;

% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% % hold on
% 
% loglog(hs,0.02*hs.^2,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^6,'--')
% hold on
%

%% Checking corrections for \partial_y (\Delta G) 

% addpath(genpath('..'))

%truev =   0.057816396686963 - 0.015079981170599i;
% truev = 1.027910069045180 + 4.450033610594533i; % 
truev = integral2(@(x,y) dens(x,y).*gradlapyonly(targ,x,y,rts,ejs),-10,10,-10,10,'AbsTol',10E-18,'RelTol',10E-18);

errs0 = hs*0;
errs1 = hs*0; 

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-10:h:10);
    src = [X(:).'; Y(:).'];
    % kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
    % val = kern{3};
    % val = val(:,:,1);
    % d1 = dens(X,Y);
    % dint = sum(val(:).*d1(:),'all');
    % errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(h,1);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
    val = kern{3};
    val = val(:,:,2);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end
K2 = errs1;


% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% % hold on
% 
% loglog(hs,0.02*hs.^2,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^6,'--')
% hold on



%% Checking corrections for \Del G

% addpath(genpath('..'))

truev = -0.361849809804212 + 1.456685781616957i;
truev = integral2(@(x,y) dens(x,y).*(hessxxonly(targ,x,y,rts,ejs)+hessyyonly(targ,x,y,rts,ejs)),-10,10,-10,10,'AbsTol',10E-18,'RelTol',10E-18);

errs0 = hs*0;
errs1 = hs*0; 

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-10:h:10);
    src = [X(:).'; Y(:).'];
    % kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
    % val = kern{2};
    % val = val(:,:,1);
    % d1 = dens(X,Y);
    % dint = sum(val(:).*d1(:),'all');
    % errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(h,1);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
    val = kern{2};
    val = val(:,:,1) + val(:,:,3);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end
K3 = errs1;


% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% hold on
% 
% loglog(hs,0.1*hs.^2,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^6,'--')
% hold on
% legend('no correction', '5 pt correction', 'h^2', 'h^6','Location','northwest') 



%% Checking corrections for G_{yy} 

% addpath(genpath('..'))

truev = -0.361849809804212 + 1.456685781616957i;
truev = integral2(@(x,y) dens(x,y).*hessyyonly(targ,x,y,rts,ejs),-10,10,-10,10,'AbsTol',10E-18,'RelTol',10E-18);

errs0 = hs*0;
errs1 = hs*0; 

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-10:h:10);
    src = [X(:).'; Y(:).'];
    % kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
    % val = kern{2};
    % val = val(:,:,1);
    % d1 = dens(X,Y);
    % dint = sum(val(:).*d1(:),'all');
    % errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(h,1);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
    val = kern{2};
    val = val(:,:,3);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end
K4 = errs1;


% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% hold on
% 
% loglog(hs,0.1*hs.^2,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^6,'--')
% hold on
% legend('no correction', '5 pt correction', 'h^2', 'h^6','Location','northwest') 

%% Checking corrections for G_{xx} 

% addpath(genpath('..'))

truev = -0.361849809804212 + 1.456685781616957i;
truev = integral2(@(x,y) dens(x,y).*hessxxonly(targ,x,y,rts,ejs),-10,10,-10,10,'AbsTol',10E-18,'RelTol',10E-18);

errs0 = hs*0;
errs1 = hs*0; 

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-10:h:10);
    src = [X(:).'; Y(:).'];
    % kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
    % val = kern{2};
    % val = val(:,:,1);
    % d1 = dens(X,Y);
    % dint = sum(val(:).*d1(:),'all');
    % errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(h,1);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
    val = kern{2};
    val = val(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end
K5 = errs1;


% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% hold on
% 
% loglog(hs,0.1*hs.^2,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^6,'--')
% hold on
% legend('no correction', '5 pt correction', 'h^2', 'h^6','Location','northwest') 


%% Checking corrections for G_{xy} 

% addpath(genpath('..'))

truev = 1.494565115435828 + 1.279414908871585i; %
truev = integral2(@(x,y) dens(x,y).*hessxyonly(targ,x,y,rts,ejs),-10,10,-10,10,'AbsTol',10E-18,'RelTol',10E-18);

errs0 = hs*0;
errs1 = hs*0; 

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-10:h:10);
    src = [X(:).'; Y(:).'];
    % kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
    % val = kern{2};
    % val = val(:,:,2);
    % d1 = dens(X,Y);
    % dint = sum(val(:).*d1(:),'all');
    % errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(h,1);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
    val = kern{2};
    val = val(:,:,2);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end
K6 = errs1;

% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% hold on
% 
% loglog(hs,0.1*hs.^4,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^6,'--')
% hold on
%
% legend('no correction', '13 pt correction', 'h^4', 'h^6','Location','northwest') 




%% Checking corrections for G 

% addpath(genpath('..'))

greenfac = @(x,y) greenvalonly(targ,x,y,rts,ejs);
truev = -3.056564167662261 - 0.163802105771794i; % 
truev = integral2(@(x,y) dens(x,y).*greenfac(x,y),-50,50,-50,50,"AbsTol",10E-18,"RelTol",10E-18);

errs0 = hs*0;
errs1 = hs*0;

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-10:h:10);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
    val = kern{1};
    d1 = dens(X,Y);
    dint = sum(val.'.*d1(:));
    errs0(ii) = abs(dint - truev) / abs(truev);

    [inds, corrs] = get_correct(h,1);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
    val = kern{1};
    dint = sum(val.*d1(:).');
    errs1(ii) = abs(dint - truev) / abs(truev);

end
K7 = errs1;

% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% hold on
% 
% loglog(hs,0.0005*hs.^4,'--')
% hold on
% 
% loglog(hs,0.0001*hs.^6,'--')
% hold on 
% 
% legend('no correction', '5 pt correction', 'h^4', 'h^6') 

%%  Checking corrections for phi 


greenfac = @(x,y) phivalonly(targ,x,y,rts,ejs);
%truev = -0.062850332948632 + 0.069122770940408i; %  integral2(@(x,y) dens(x,y).*greenfac(x,y),-50,50,-50,50,"AbsTol",0,"RelTol",10E-14);
truev =  -1.180843688712207 - 0.251526644367094i; % 
truev = integral2(@(x,y) dens(x,y).*greenfac(x,y),-30,30,-30,30,"AbsTol",10E-18,"RelTol",10E-18);

errs0 = hs*0;
errs1 = hs*0;

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-10:h:10);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
    val = kern{4};
    d1 = dens(X,Y);
    dint = sum(val.'.*d1(:));
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(h,1);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
    val = kern{4};
    dint = sum(val.'.*d1(:));
    errs1(ii) = abs(dint - truev) / abs(truev);

end
K8 = errs1;

% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% hold on
% 
% loglog(hs,0.0005*hs.^5,'--')
% hold on
% 
% loglog(hs,0.0001*hs.^6,'--')
% hold on 
% 
% legend('no correction', '5 pt correction', 'h^5', 'h^6') 
%%

f = figure(1); clf
f.Position = [1 1 634 386];

n = 8;
C = orderedcolors("gem");            % or get(groot,'defaultAxesColorOrder')
C(8,:) = [0.5111 0.2000 0.1778];     % brick
colororder(C)
mk = {'o','s','^','d','v','x','+','p'};   % 'Marker', mk{k}

hs = [2 1 0.5 0.2 0.1 0.05 0.025 0.025/2 0.025/4 0.025/8];

loglog(hs,K1,'o-','DisplayName','$K_1$','LineWidth',0.8,'Marker', mk{1})
hold on

loglog(hs,K2,'o-','DisplayName','$K_2$','LineWidth',0.8,'Marker', mk{2})
hold on

loglog(hs,K3,'o-','DisplayName','$K_3$','LineWidth',0.8,'Marker', mk{3})
hold on

loglog(hs,K4,'o-','DisplayName','$K_4$','LineWidth',0.8,'Marker', mk{4})
hold on

loglog(hs,K5,'o-','DisplayName','$K_5$','LineWidth',0.8,'Marker', mk{1})
hold on

loglog(hs,K6,'o-','DisplayName','$K_6$','LineWidth',0.8,'Marker', mk{2})
hold on

loglog(hs,K7,'o-','DisplayName','$K_7$','LineWidth',0.8,'Marker', mk{3})
hold on

loglog(hs,K8,'o-','DisplayName','$K_8$','LineWidth',0.8,'Marker', mk{4})
hold on

hs = [0.5 0.2 0.1 0.05 0.025 0.025/2 0.025/4];

loglog(hs,hs.^6,'k--','LineWidth',0.8,'DisplayName','$h^6$')

legend('Interpreter','latex','Location','eastoutside')

ylabel('Relative error')
xlabel('$h$','Interpreter','latex')

xlim([2E-3 3])
ylim([1E-16 2])

set(gca, 'FontSize',12)
fontname(gcf, 'CMU Serif')


%%

saveas(gcf,'intconv.fig','fig')
exportgraphics(gcf,'intconv.pdf','ContentType','vector')

