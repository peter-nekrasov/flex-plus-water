function [coefs,H] = pit(X,Y,amp,width,w)

    E = 7*10^9;
    nu = 0.33;
    H0 = 1;
    rhow = 1025;
    rhoi = 917;
    g = 9.8;

    a0 = E*H0^3/(12*(1-nu^2));
    b0 = (rhoi*H0*w^2 - rhow*g);
    g0 = -w^2*rhow;
    gbar = 0;
    Hbar = amp*exp(-(X.^2 + Y.^2)/(2*width^2));
    
    s = 0.0002;
    H = H0+amp/2+amp/2*erf(s*(-X.^2 - Y.^2 + width^2));


    alpha = E*H.^3/(12*(1-nu^2));
    beta = (rhoi*H*w^2 - rhow*g);

    abar = alpha - a0;
    bbar = beta - b0;


    Hx = -2*amp*s*X.*exp(-s^2*(-X.^2 - Y.^2 + width^2).^2)/sqrt(pi);
    Hy = -2*amp*s*Y.*exp(-s^2*(-X.^2 - Y.^2 + width^2).^2)/sqrt(pi);

    Hxx = -2*amp*s*exp(-s^2*(-X.^2 - Y.^2 + width^2).^2)/sqrt(pi) - 8*amp*s^3*X^2.*(-X.^2 - Y.^2 + width^2).*exp(-s^2*(-X.^2 - Y.^2 + width^2).^2)/sqrt(pi);
    Hxy = -8*amp*s^3*X.*Y.*(-X.^2 - Y.^2 + width^2).*exp(-s^2*(-X.^2 - Y.^2 + width^2).^2)/sqrt(pi);
    Hyy = -2*amp*s*exp(-s^2*(-X.^2 - Y.^2 + width^2).^2)/sqrt(pi) - 8*amp*s^3*Y^2.*(-X.^2 - Y.^2 + width^2).*exp(-s^2*(-X.^2 - Y.^2 + width^2).^2)/sqrt(pi);

    ax = 3*E*H.^2.*Hx/(12*(1-nu^2));
    ay = 3*E*H.^2.*Hy/(12*(1-nu^2));

    axx = 6*E*H.*Hx.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hxx/(12*(1-nu^2));
    axy = 6*E*H.*Hx.*Hy/(12*(1-nu^2)) + 3*E*H.^2.*Hxy/(12*(1-nu^2));
    ayy = 6*E*H.*Hy.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hyy/(12*(1-nu^2));

    coefs = {a0,abar,b0,bbar,g0,gbar,ax,ay,axx,axy,ayy,nu}; 

end