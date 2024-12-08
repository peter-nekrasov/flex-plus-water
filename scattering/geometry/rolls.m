function [coefs,H] = rolls(X,Y,xl,xr,yl,yr,amp,width,freq)

    E = 7*10^9;
    nu = 0.33;
    rhow = 1025;
    rhoi = 917;
    g = 9.8;
    w = freq; % frequency

    freb = 0.35;

    H0 = -rhow/(rhoi-rhow)*freb ;

    [psi,psix,psiy,psixx,psixy,psiyy]  = bumpfunc(X,Y,xl,xr,yl,yr);

    S = (amp/2+amp/2*sin(2*pi*X/width)).*psi + freb ;
    B = rhoi/(rhoi-rhow)*S;
    H = S - B;

    Sx = amp*pi/width*cos(2*pi*X/width).*psi + (amp/2+amp/2*sin(2*pi*X/width)).*psix;
    Sy = (amp/2+amp/2*sin(2*pi*X/width)).*psiy;

    Sxx = -2*amp*pi^2/width^2*sin(2*pi*X/width).*psi + ...
        amp*pi/width*cos(2*pi*X/width).*psix + ...
        amp*pi/width*cos(2*pi*X/width).*psix + ...
        (amp/2+amp/2*sin(2*pi*X/width)).*psixx;
    Sxy = amp*pi/width*cos(2*pi*X/width).*psiy + (amp/2+amp/2*sin(2*pi*X/width)).*psixy;
    Syy = (amp/2+amp/2*sin(2*pi*X/width)).*psiyy;

    Hx = rhow/(rhow - rhoi).*Sx;
    Hy = rhow/(rhow - rhoi).*Sy;
    Hxx = rhow/(rhow - rhoi).*Sxx;
    Hxy = rhow/(rhow - rhoi).*Sxy;
    Hyy = rhow/(rhow - rhoi).*Syy;

    a0 = E*H0^3/(12*(1-nu^2));
    b0 = (rhoi*H0*w^2 - rhow*g);
    g0 = -w^2*rhow;
    gbar = 0;

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