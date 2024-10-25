function out = gaussian(X,Y,H0,width)

    E = 7E9; % Young's modulus (Pa)
    nu = 0.3; % Poisson ratio

    H = H0 + exp(-(X.^2 + Y.^2)/(2*width^2));

    D = E*H^3/(12*(1-nu^2));

    Dx = - 3*X.*D/width^2;
    Dy = - 3*Y.*D/width^2;

    Dxx = - 3*D/width^2 - 3*X.*Dx/width^2;
    Dxy = - 3*X.*Dy/width^2;
    Dyy = - 3*D/width^2 - 3*Y.*Dy/width^2;

    out = {H,D,Dx,Dy,Dxx,Dxy,Dyy};

end