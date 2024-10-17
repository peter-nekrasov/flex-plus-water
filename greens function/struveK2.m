function [val,grad] = struveK2(rhoj,src,targ)
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
    
    ilow = (imag(rhoj) < 0);

    if ilow
        rhoj = conj(rhoj); % flip rhoj up into upper half plane
    end

    zt = r*rhoj;
    [cr0,cr1] = struveR(zt);

    [h0,gradh0] = helmdiffgreen(rhoj,src,targ);
    h0(zt == 0) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));
    h0 = -4i*h0;
    gradh0 = -4i*gradh0;
    gradh0x = gradh0(:,:,1);
    gradh0y = gradh0(:,:,2);


    val = -1i*cr0+1i*h0;
    gradx = 1i*rhoj*cr1.*xt./r+1i*gradh0x;
    grady = 1i*rhoj*cr1.*yt./r+1i*gradh0y;

    gradx(r == 0) = 0;
    grady(r == 0) = 0;

    if ilow
        val = conj(val);
        gradx = conj(gradx);
        grady = conj(grady);
    end

    grad = cat(3,gradx,grady);

end