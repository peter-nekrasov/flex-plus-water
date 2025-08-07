function [proxy_dict] = init_proxy_ann_dict(opts)
    % if ~isfield(opts,'k'), opts.k  = 8; end
    % if ~isfield(opts,'krad'), opts.krad  = 2; end
    % if ~isfield(opts,'rat'), opts.rat  = 1.5; end
    if ~isfield(opts,'k'), opts.k  = 10; end
    if ~isfield(opts,'krad'), opts.krad  = 10; end
    if ~isfield(opts,'rat'), opts.rat  = 1.5; end
    % if ~isfield(opts,'krad'), opts.krad  = 40; end
    % if ~isfield(opts,'rat'), opts.rat  = 1.5; end
    if ~isfield(opts,'dim'), opts.dim = 3; end


    RR = 5/2;
    if opts.dim == 3
        [ppts,wpts] = prx_points_cube(RR,RR*opts.rat,opts.k,opts.krad);
    else
        [ppts,wpts] = prx_points_square(RR,RR*opts.rat,opts.k,opts.krad);
    end
    ppts = ppts-0.5;

    proxy_dict = [];
    proxy_dict.proxy = ppts;
    proxy_dict.weigt = wpts;
    proxy_dict.norms = NaN*ppts;
end


function [xprox,wprox] = prx_points_cube(Lprxin,Lprxout,k,krad)

[x_leg,w_leg] = legpts(k, [0,1]);

[x_rad,w_rad] = legpts(krad, [Lprxin,Lprxout]);

xface = [repmat(x_leg,k,1),reshape(repmat(x_leg',k,1),[],1)]';
wface = repmat(w_leg,1,k).*(reshape(repmat(w_leg,k,1),[],1)');
xface3 = [xface;zeros(1,k*k)];

xcube = [xface3,xface3([1,3,2],:),xface3([3,1,2],:)];
xcube = [xcube,1-xcube];
wcube = repmat(wface,1,6);

xprox=arrayfun(@(l) (xcube-0.5)*l+0.5,x_rad,'Uniform',0);
xprox = cat(2,xprox{:});
wprox=arrayfun(@(l,w) .5*wcube*l*l*w,x_rad,w_rad','Uniform',0);
wprox = cat(2,wprox{:});
% nprox = numel(wprox);

end

function [xprox,wprox] = prx_points_square(Lprxin,Lprxout,k,krad)

[x_leg,w_leg] = legpts(k, [0,1]);

[x_rad,w_rad] = legpts(krad, [Lprxin,Lprxout]);

xface = x_leg(:)';
wface = w_leg;
xface2 = [xface;zeros(1,k)];

xsquare = [xface2, xface2([2,1],:)];
xsquare = [xsquare,1-xsquare];
wsquare = repmat(wface,1,4);

xprox=arrayfun(@(l) (xsquare-0.5)*l+0.5,x_rad,'Uniform',0);
xprox = cat(2,xprox{:});
wprox=arrayfun(@(l,w) .5*wsquare*l*l*w,x_rad,w_rad','Uniform',0);
wprox = cat(2,wprox{:});
% nprox = numel(wprox);

end
