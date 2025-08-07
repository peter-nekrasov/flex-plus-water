function kerns = gen_fft_kerns_helm(kerns,sz,ind)
        
    GS = reshape(kerns{1},sz);

    [zi,zj] = ind2sub(sz,ind);

    GS = circshift(GS,[zi,zj]);
    kerns = fft2(GS);


end