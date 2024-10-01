%test001
%
% test of root finding for the polynomial 
%
% f = z^5 + d1 z + d0 
%
% which is in Bring-Jerrard form
%
% We are also interested in the partial fraction expansion
% of g = z/f 
% 

beta = 1.1 + 1i*0.0001;
gamma = -(1.1)^(5/4)*((1/5)^(1/4)-(1/5)^(5/4));
d1 = -(beta);
d0 = gamma;

f = @(z) z.^5 + d1*z + d0;
fp = @(z) 5*z.^4 + d1;

C = bring_companion(d1,d0);

rts1 = eig(C);
rts2 = bring_refine_rts(rts1,d1,d0);

assert(max(abs(f(rts2)./fp(rts2))) < 30*eps(1))

rts2

%%
ejs = 1./fp(rts2); % partial fraction coefficients for 1/f

% properties of ejs and roots via elementary algebra and residue calculus 
abs(sum(rts2))
abs(prod(rts2)/d0+1)
abs(sum(ejs.*rts2))
abs(sum(ejs.*rts2.^2))
abs(sum(ejs.*rts2.^3))
abs(sum(ejs.*rts2.^4)-1)

zrand = randn(1000,1,'like',1.0+1.0*1i);

g = @(z) z./f(z);
pfac = @(z) sum(ejs(:).'.*rts2(:).'./(z(:)-rts2(:).'),2);

max(abs( g(zrand) - pfac(zrand)))/max(abs(g(zrand)))

%% conditioning of bring_companion matrix/count of real roots

d0s = linspace(-10,10,101);
d1s = linspace(-10,10,101);

conds = zeros(length(d1s),length(d0s));
fpmins = zeros(length(d1s),length(d0s));
realcount = zeros(length(d1s),length(d0s));

for j = 1:length(d0s)
    for i = 1:length(d1s)
        d0 = d0s(i);
        d1 = d1s(j);
        C = bring_companion(d1,d0);
    
        conds(i,j) = cond(C);

        f = @(z) z.^5 + d1*z + d0;
        fp = @(z) 5*z.^4 + d1;
        rts1 = eig(C);
        rts2 = bring_refine_rts(rts1,d1,d0);
        realcount(i,j) = nnz(abs(imag(rts2))< 1e-12);

        fpmins(i,j) = min(abs(fp(rts2)));
        
    end
end

figure(1)
clf
h = imagesc(d0s,d1s,realcount);
colorbar

%% critical gamma as a function of beta
% d1 = -beta, d0 = gamma 

betas = linspace(0,10,100);

gammalow = 0; gammahigh = 20;
gammas = zeros(size(betas));

ntry = 54;
for j = 2:length(betas)
    a = gammalow; b = gammahigh;
    beta = betas(j);
    for ii = 1:ntry
        c = (a+b)/2;
        d0 = c;
        d1 = -beta;
        C = bring_companion(d1,d0);
        rts1 = eig(C);
        rts2 = bring_refine_rts(rts1,d1,d0);
        realcount = nnz(and(abs(imag(rts2))< 1e-12, real(rts2)>0))
        if realcount == 1
            disp('interesting')
        end
        if realcount == 2
            a = c;
        else
            b = c;
        end
    end
    gammas(j) = a;
end

%%

gammaformula = betas.^(5/4)*(1/5^(1/4)-1/5^(5/4));
clf
plot(betas,gammas,'r')
hold on
plot(betas,gammaformula,'g-.')


function C = bring_companion(d1,d0)
%
% returns the companion matrix of z^5 + d1 z + d0
%

C = diag(ones(4,1),-1);
C(:,end) = [-d0; -d1; 0; 0; 0];

end

function rts2 = bring_refine_rts(rts1,d1,d0,nstep)
%
% applies nstep steps of Newton to each root
%

if nargin < 4
    nstep = 3;
end

f = @(z) z.^5 + d1*z + d0;
fp = @(z) 5*z.^4 + d1;

rts2 = rts1;

for j = 1:nstep
    rts2 = rts2 - f(rts2)./fp(rts2);
end
end