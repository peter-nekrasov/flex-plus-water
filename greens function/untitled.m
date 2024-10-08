
beta = 1.1;
gamma = (1.1)^(5/4)*((1/5)^(1/4)-(1/5)^(5/4));

beta = 1.1+1e-2*1i;
gamma = 0.2;
epsilon = 0.98e-2;

f = @(z) z.^5 - beta*z + gamma;
fp = @(z) 5*z.^4 - beta;

f2 = @(z) z.^7 - beta*z.^3 + gamma*(z.^2 - epsilon^2*0.5);
fp2 = @(z) 7*z.^6 - 3*beta*z.^2 + 2*gamma*z;

pp = [1 0 0 0 -beta gamma 0 -gamma*0.5*epsilon^2];
pp0 = [1 0 0 0 -beta gamma];

rts = roots(pp);
rts0 = roots(pp0);

abs(f2(rts))
abs(f(rts0))

rts
rts0