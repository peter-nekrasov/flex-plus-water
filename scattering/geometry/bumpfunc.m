function [psi,psix,psiy,psixx,psixy,psiyy] = bumpfunc(X,Y,xl,xr,yl,yr)

    s = 0.008;

    psi =  ( erf(s*(X - xl)) + erf(s*(-X + xr))).*( erf(s*(Y - yl)) + erf(s*(-Y + yr)))/4;

    psix = 1/2/sqrt(pi)*s*( exp(-s^2*(X - xl).^2) - exp(-s^2*(X - xr).^2)).*( erf(s*(Y - yl)) + erf(s*(-Y + yr)));
    psiy = 1/2/sqrt(pi)*s*( exp(-s^2*(Y - yl).^2) - exp(-s^2*(Y - yr).^2)).*( erf(s*(X - xl)) + erf(s*(-X + xr)));
    
    psixx = 1/sqrt(pi)*s^3*((-X + xl).*exp(-(s*(X - xl)).^2) + (X - xr).*exp(-(s*(-X + xr)).^2)).*( erf(s*(Y - yl)) + erf(s*(-Y + yr)));
    psiyy = 1/sqrt(pi)*s^3*((-Y + yl).*exp(-(s*(Y - yl)).^2) + (Y - yr).*exp(-(s*(-Y + yr)).^2)).*( erf(s*(X - xl)) + erf(s*(-X + xr)));
    
    psixy = s^2/pi*( exp(-(s*(X - xl)).^2) - exp(-(s*(-X + xr)).^2)).*( exp(-(s*(Y - yl)).^2) - exp(-(s*(-Y + yr)).^2)); 

end