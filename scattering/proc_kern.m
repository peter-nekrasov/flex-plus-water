function kern_struct = proc_kern(kern_struct,h)

    val = kern_struct{1};
    hessxx = kern_struct{2};
    hessxy = kern_struct{3};
    hessyy = kern_struct{4};
    gradlapx = kern_struct{5};
    gradlapy = kern_struct{6};
    phi = kern_struct{7};

    %%% PERFORM DIAGONAL CORRECTIONS HERE
    % NEED TO KNOW CONSTANT PART 

    
    for ii = 1:numel(kern_struct)
        kern = kern_struct{ii};
        kernhat = [kern, flip(kern(1:end,2:end),2); ...
        flip(kern(2:end,1:end)), flip(flip(kern(2:end,2:end)),2)];
        kern_struct{ii} = fft2(kernhat)*h*h;
    end


end