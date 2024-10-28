function M = eval_kern(src,targ,func,inds,corrs)
% Evaluates kernels and adds corrections to the source part based on inds

    kerns = func(src,targ);

    if (numel(kerns) ~= numel(inds)) | (numel(inds) ~= numel(corrs))
        error('output of func must be the same length as inds and corrs')
    end

    M = cell(1,numel(kerns));

    for ii = 1:numel(kerns)
        kern = kerns(ii);
        kern(inds) = kern(inds) + corrs;
        M{ii} = kern;
    end

end