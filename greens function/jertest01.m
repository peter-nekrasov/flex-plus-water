beta = 2;
gamma = - 1 - 0.2i;
[rts,ejs] = find_roots(beta,gamma);

rs = [120; 150; 200];
vals = rs*0;

% rhos = 0:0.1:1500;
% vals = integrand(rhos,pi/2);
% 
% loglog(rhos,abs(vals))
% hold on
% loglog(rhos(300:end),0.01*rhos(300:end).^(-1/2))
% xlabel('\rho')
% legend('abs','r^{-1/2}')


for ii = 1:numel(rs)
    rtargy = rs(ii);
    rtargx = 0;
    integrand = @(rho,theta) greenvalonly([5;5],rtargx-rho.*cos(theta),rtargy - rho.*sin(theta),rts,ejs)/(2*pi);
    vals(ii) = integral2(@(rho,theta) integrand(rho,theta), 0,300+rtargy,0,2*pi)
end


%% calculate finite difference term

h = 0.001;
rs = [    10
    20
    40
    60
    80
   100
   120
   150];

% finite difference stencils 
% d4 - partial_{xxxx} (6th order) % improve order?
d4 = zeros(9,1);
d4(1) = 7/240;	
d4(2) = -2/5;
d4(3) = 169/60;
d4(4) = -122/15;
d4(5) = 91/8;
d4(6) = -122/15;
d4(7) = 169/60;
d4(8) = -2/5;
d4(9) = 7/240;

% d2 - partial_{xx} (8th order)
d2 = zeros(9, 1);
d2(1) = -1/560;
d2(2) = 8/315;
d2(3) = -1/5;
d2(4) = 8/5;
d2(5) = -205/72;
d2(6) = 8/5;
d2(7) = -1/5;
d2(8) = 8/315;
d2(9) = -1/560;

bilap = zeros(9);
bilap(:,5) = d4;
bilap(5,:) = bilap(5,:) + d4.';
bilap = bilap + 2*(d2*d2.');
bilap = bilap / h^4;

src = [5;5];

firstterm = rs*0;

for ii = 1:numel(rs)
    rtemp = rs(ii);
    [X,Y] = meshgrid(-4*h:h:4*h, (rtemp-4*h):h:(rtemp+4*h));
    targ = [X(:).'; Y(:).'];
    gs = green(src,targ,rts,ejs);
    gval = gs{1};
    gval = reshape(gval, size(X));
    firstterm(ii) = (sum(bilap.*gval,'all') - beta*gval(5:5))./gamma;
end

firstterm



%%
t = tiledlayout(1,2);

rs2 = sqrt((rs - 5).^2 + 5.^2);


nexttile
plot(rs2,(abs(vals)),'o-')
ylabel('')
xlabel('R')

nexttile
semilogy(rs2,(abs(vals)),'o-')
xlabel('R')

title(t,'Integral of G_S against 1/|r| for complex \gamma')
return


function val = greenvalonly(targ,x,y,rts,ejs)

    src = [x(:).'; y(:).'];
    gf = green(src,targ,rts,ejs);
    val = reshape(gf{1},size(x));
    

end