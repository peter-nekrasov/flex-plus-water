beta = 1.1;
gamma = (1.1)^(5/4)*((1/5)^(1/4)-(1/5)^(5/4));

beta = 1.1;
gamma = 0.2;
epsilon = 0.98e-2;

f = @(z) z.^5 - beta*z + gamma;
fp = @(z) 5*z.^4 - beta;

f2 = @(z) epsilon*z.^7 + f(z);
fp2 = @(z)epsilon*7*z.*6 + fp(z);

pp = [1i*epsilon 0 1 0 1i*epsilon 0 -beta gamma];
pp0 = [1 0 0 0 -beta gamma];

rts = roots(pp);
rts0 = roots(pp0);

abs(f2(rts))
abs(f(rts0))

rts
rts0
