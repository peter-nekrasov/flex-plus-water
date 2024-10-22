%% 1d trapezoid rule integration test
% for a smooth function, convergence is super-algebraic

hs = 1./2.^(-10:15);
dens = @(x) x.^2; % exp(-(x.^2)/(2*5^2)); % test density 
errs = hs*0;

trueint = 2*50^3/3; % 12.53314137315500251207863542266483733193

for i = 1:numel(hs)
    h = hs(i);
    xs = -50:h:50;
    errs(i) = abs(trueint - sum(h*dens(xs)))/trueint;
end

loglog(hs,errs)
hold on

loglog(hs,hs.^1)

