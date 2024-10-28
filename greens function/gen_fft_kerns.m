function kerns = gen_fft_kerns(kerns,sz,ind)
        
    kerns{1} = reshape(kerns{1},sz);
    kerns{4} = reshape(kerns{4},sz);

    [zi,zj] = ind2sub(sz,ind);

    kernhat = circshift(kerns{1},[zi,zj]);
    kerns{1} = fft2(kernhat);

    kernhat = circshift(kerns{4},[zi,zj]);
    kerns{4} = fft2(kernhat);

end