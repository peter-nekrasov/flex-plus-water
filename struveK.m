function [ck0,ck1] = struveK(rhoj, z)
    % where rhoj is a complex number and z is an array of real numbers

    sz = size(z);
    z = z(:).';
    
    ilow = (imag(rhoj) < 0);

    if ilow
        rhoj = conj(rhoj); % flip rhoj up into upper half plane
    end

    zt = z*rhoj;

    [cr0,cr1] = struveR(zt);

    src = [0; 0]; targ = [ zt / rhoj ; 0*z];
    [h0,~] = helmdiffgreen(rhoj,src,targ);
    h0 = -4i*h0.' ;
    h1 = besselh(1,zt); % need to use evaluator here for H1

    ck0 = -1i*cr0+1i*h0;
    ck1 = -1i*cr1+1i*h1;

    if ilow
        ck0 = conj(ck0);
        ck1 = conj(ck1);
    end

    ck0 = reshape(ck0,sz);
    ck1 = reshape(ck1,sz);

end
