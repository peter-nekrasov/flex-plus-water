function [s,sp,spp] = sigmoid(X)
    a = 0.01;
    s = 1./(1+exp(-a*X));
    sp = a*exp(a*X)./(exp(a*X)+1).^2;
    spp = -a^2*(exp(a*X)-1).*exp(a*X)./(exp(a*X)+1).^3;
end