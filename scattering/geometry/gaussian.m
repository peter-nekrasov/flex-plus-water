function out = gaussian(X,Y)

    width = 5;

    val = exp(-(X.^2 + Y.^2)/(2*width^2));

    gradx = - X.*val/width^2;
    grady = - Y.*val/width^2;

    hessxx = - val/width^2 - X.*gradx/width^2;
    hessxy = - X.*grady/width^2;
    hessyy = - val/width^2 - Y.*grady/width^2;

    out = {val,gradx,grady,hessxx,hessxy,hessyy};

end