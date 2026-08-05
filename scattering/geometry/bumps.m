function [coefs,H] = bumps(X,Y,xmin,xmax,amp,width,w)

    rng(45);

    E = 7*10^9;
    nu = 0.33;
    H0 = 5;
    rhow = 1025;
    rhoi = 917;
    g = 9.8;

    a0 = E*H0^3/(12*(1-nu^2));
    b0 = (rhoi*H0*w^2 - rhow*g);
    g0 = -w^2*rhow;
    gbar = 0;

    [Xg,Yg] = meshgrid(xmin:width:xmax);

    Xg = Xg(:);
    Yg = Yg(:);

    H = H0;
    Hx = 0;
    Hy = 0;
    Hxx = 0;
    Hxy = 0;
    Hyy = 0;

    for ii = 1:numel(Xg)

        Xo = Xg(ii);
        Yo = Yg(ii);

        Hbar = (2*rand-1)*amp*exp(-((X-Xo).^2 + (Y-Yo).^2)/(2*width^2));

        H = H + Hbar;

        tpHx = -(X-Xo).*Hbar/width^2;
        tpHy = -(Y-Yo).*Hbar/width^2;
    
        tpHxx = - Hbar/width^2 - (X-Xo).*tpHx/width^2;
        tpHxy = - (X-Xo).*tpHy/width^2;
        tpHyy = - Hbar/width^2 - (Y-Yo).*tpHy/width^2;

        Hx = Hx + tpHx;
        Hy = Hy + tpHy;

        Hxx = Hxx + tpHxx;
        Hxy = Hxy + tpHxy;
        Hyy = Hyy + tpHyy;

    end

    alpha = E*H.^3/(12*(1-nu^2));
    beta = (rhoi*H*w^2 - rhow*g);

    abar = alpha - a0;
    bbar = beta - b0;

    ax = 3*E*H.^2.*Hx/(12*(1-nu^2));
    ay = 3*E*H.^2.*Hy/(12*(1-nu^2));

    axx = 6*E*H.*Hx.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hxx/(12*(1-nu^2));
    axy = 6*E*H.*Hx.*Hy/(12*(1-nu^2)) + 3*E*H.^2.*Hxy/(12*(1-nu^2));
    ayy = 6*E*H.*Hy.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hyy/(12*(1-nu^2));

    coefs = {a0/E,abar/E,b0/E,bbar/E,g0/E,gbar/E,ax/E,ay/E,axx/E,axy/E,ayy/E,nu}; 

end