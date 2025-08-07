function [Kpxy,lst] = proxyfun_kern_rskelf(x,slf,lst,proxy_dict,l,ctr,kernout,kernin,srcinfo)
    l = max(l);
    proxy = proxy_dict.proxy;
    weigt = l*proxy_dict.weigt;
    norms = proxy_dict.norms;
    % Shift and scale precomputed proxy surface
    % pxy = bsxfun(@plus, proxy*l, ctr');
    pxy = proxy*l + ctr;
    % set up sourceinfo and targetinfo
    % srcinfo = [];
    % srcinfo.r = x(:,slf);
    % srcinfo.n = rn(:,slf);
    % w_sqrt = sqrt(wts(slf).');

    srcuse = [];
    srcuse.r = x(:,slf);
    if isfield(srcinfo,'n')
        srcuse.n = srcinfo.n(:,slf);
    end
    if isfield(srcinfo,'d')
        srcuse.d = srcinfo.d(:,slf);
    end
    if isfield(srcinfo,'d2')
        srcuse.d2 = srcinfo.d2(:,slf);
    end

    w_sqrt = (srcinfo.wts(slf).');
    
    
    
    pxyzinfo = [];
    pxyzinfo.r = pxy;
    pxyzinfo.n = norms;
    
    % Kpxy1 = bsxfun(@times,belpde.kern(zk,pxyzinfo, ...
    %   srcinfo,kernin).',w_sqrt);
    Kpxy1 = bsxfun(@times,kernin(pxyzinfo,srcuse).',w_sqrt);

    Kpxy1 = bsxfun(@times, weigt.', Kpxy1);
    % Kpxy2 = bsxfun(@times,belpde.kern(zk,srcinfo, ...
    %   pxyzinfo,kernout),w_sqrt);
    Kpxy2 = bsxfun(@times,kernout(srcuse,pxyzinfo),w_sqrt);
    Kpxy2 = bsxfun(@times, weigt.', Kpxy2);
    
    Kpxy = [Kpxy1;Kpxy2];%;ones(1,numel(slf))];
    dxyz = abs(x(1:2,lst)-ctr(1:2))/l;
    lst = lst(prod(dxyz < 2.5,1)>1-1e-1); %Cubical proxy
    % lst = lst(prod(dxyz < 4,1)>1-1e-1); %Cubical proxy

end
