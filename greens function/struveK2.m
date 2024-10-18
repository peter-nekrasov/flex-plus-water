function [val,grad,hess] = struveK2(rhoj,src,targ)
    % where rhoj is a complex number and z is an array of real numbers
    % the hankel part of this function has log subtracted

    [~,ns] = size(src);
    [~,nt] = size(targ);

    xs = repmat(src(1,:),nt,1);
    ys = repmat(src(2,:),nt,1);

    xt = repmat(targ(1,:).',1,ns);
    yt = repmat(targ(2,:).',1,ns);

    dx = xt-xs;
    dy = yt-ys;
    
    dx2 = dx.*dx;
    dy2 = dy.*dy;

    r2 = dx2 + dy2;
    r = sqrt(r2);
    r3 = r.^3;
    
    ilow = (imag(rhoj) < 0);

    if ilow
        rhoj = conj(rhoj); % flip rhoj up into upper half plane
    end

    zt = r*rhoj;
    [cr0,cr1] = struveR(zt);

    [h0,gradh0,hessh0] = helmdiffgreen(rhoj,src,targ);
    
    h0x = gradh0(:,:,1);
    h0y = gradh0(:,:,2);

    h0xx = hessh0(:,:,1);
    h0xy = hessh0(:,:,2);
    h0yy = hessh0(:,:,3);

    h0(zt == 0) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));

    h0 = -4i*h0;

    h0x = -4i*h0x;
    h0y = -4i*h0y;

    h0xx = -4i*h0xx;
    h0xy = -4i*h0xy;
    h0yy = -4i*h0yy;

    val = -1i*cr0+1i*h0;
    gradx = 1i*rhoj*cr1.*xt./r+1i*h0x;
    grady = 1i*rhoj*cr1.*yt./r+1i*h0y;
    hessxx = 1i*rhoj^2*xt.^2./r2.*cr0-1i*rhoj*(xt.^2 - yt.^2)./r3.*cr1+1i*h0xx; 
    hessxy = 1i*rhoj^2*xt.*yt./r2.*cr0-2i*rhoj*xt.*yt./r3.*cr1+1i*h0xy; 
    hessyy = 1i*rhoj^2*yt.^2./r2.*cr0+1i*rhoj*(xt.^2 - yt.^2)./r3.*cr1+1i*h0yy; 

    gradx(r == 0) = 0;
    grady(r == 0) = 0;

    hessxx(r == 0) = 1i*rhoj^2*cr0(r == 0) -4i*1i*rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-1.5-log(2));
    hessxy(r == 0) = 0;
    hessyy(r == 0) = 1i*rhoj^2*cr0(r == 0) -4i*1i*rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-1.5-log(2));


    % add corrections to second derivatives

    if ilow
        val = conj(val);
        gradx = conj(gradx);
        grady = conj(grady);
        hessxx = conj(hessxx);
        hessxy = conj(hessxy);
        hessyy = conj(hessyy);
    end

    grad = cat(3,gradx,grady);
    hess = cat(3,hessxx,hessxy,hessyy);

end