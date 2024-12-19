function [coefs,H] = dipole(X,Y,amp,width,w)

    E = 7*10^9;
    nu = 0.33;
    H0 = 20;
    rhow = 1025;
    rhoi = 917;
    g = 9.8;

    a0 = E*H0^3/(12*(1-nu^2));
    b0 = (rhoi*H0*w^2 - rhow*g);
    g0 = -w^2*rhow;
    gbar = 0;

    Hbar = amp*X.*exp(-(X.^2 + Y.^2)/(2*width^2));
    H = H0 + Hbar;

    alpha = E*H.^3/(12*(1-nu^2));
    beta = (rhoi*H*w^2 - rhow*g);

    abar = alpha - a0;
    bbar = beta - b0;

    Hx = amp*exp(-(X.^2 + Y.^2)/(2*width^2)) - X.*Hbar/width^2;
    Hy = -Y.*Hbar/width^2;
  
    Hxx = -2*Hbar/(width^2) - X.*Hx/width^2;
    Hxy = -Y.*Hx/width^2;
    Hyy = -Hbar/width^2 - Y.*Hy/width^2;

    ax = 3*E*H.^2.*Hx/(12*(1-nu^2));
    ay = 3*E*H.^2.*Hy/(12*(1-nu^2));

    axx = 6*E*H.*Hx.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hxx/(12*(1-nu^2));
    axy = 6*E*H.*Hx.*Hy/(12*(1-nu^2)) + 3*E*H.^2.*Hxy/(12*(1-nu^2));
    ayy = 6*E*H.*Hy.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hyy/(12*(1-nu^2));

    coefs = {a0/E,abar/E,b0/E,bbar/E,g0/E,gbar/E,ax/E,ay/E,axx/E,axy/E,ayy/E,nu}; 

end