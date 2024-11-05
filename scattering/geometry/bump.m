function [coefs] = bump(X,Y,amp,width)

    a0 = 3; 
    b0 = 5; 
    g0 = -1; 
    nu = 0.33;

    abar = amp*exp(-(X.^2 + Y.^2)/(2*width^2));
    bbar = amp*exp(-(X.^2 + Y.^2)/(2*width^2));
    gbar = 0;

    ax = - X.*abar/width^2;
    ay = - Y.*abar/width^2;

    axx = - abar/width^2 - X.*ax/width^2;
    axy = - X.*ay/width^2;
    ayy = - abar/width^2 - Y.*ay/width^2;

    coefs = {a0,abar,b0,bbar,g0,gbar,ax,ay,axx,axy,ayy,nu}; 

end