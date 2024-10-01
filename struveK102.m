% function [ck0,ck1] = struveK102(z)
% 
%     sz = size(z);
%     z = z(:);
% 
%     [ck0,ck1] = struve102(-z);
% 
%     h0 = besselh(0,-z);
%     h1 = besselh(1,-z);
% 
%     ck0 = -1i*ck0+1i*h0;
%     ck1 = -1i*ck1+1i*h1; % needs to be checked
% 
%     ck0 = reshape(ck0,sz);
%     ck1 = reshape(ck1,sz);
% 
% end


function [ck0,ck1] = struveK102(z)

    sz = size(z);
    z = -z(:);

    iupp = find(imag(z) >= 0);
    ilow = find(imag(z) <= 0);

    zt = z;
    zt(ilow) = -z(ilow);

    ck0 = zeros(size(z));
    ck1 = zeros(size(z));

    [ck0(iupp),ck1(iupp)] = struve102(zt(iupp));
    [ck0(ilow),ck1(ilow)] = struve102(zt(ilow));

    h0 = besselh(0,zt);
    h1 = besselh(1,zt);

    ck0(iupp) = -1i*ck0(iupp)+1i*h0(iupp);
    ck1(iupp) = -1i*ck1(iupp)+1i*h1(iupp);

    ck0(ilow) =  1i*ck0(ilow)+1i*h0(ilow);
    ck1(ilow) = -1i*ck0(ilow)-1i*h1(ilow);

    ck0 = reshape(ck0,sz);
    ck1 = reshape(ck1,sz);

end
