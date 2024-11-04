function kerns = gen_fft_kerns(kerns,sz,ind)
        
    GS = reshape(kerns{1},sz);
    Gphi = reshape(kerns{4},sz);

    hessGS = kerns{2};
    gradlapGS = kerns{3};

    [zi,zj] = ind2sub(sz,ind);

    GS = circshift(GS,[zi,zj]);
    kerns{1} = fft2(GS);

    Gphi = circshift(Gphi,[zi,zj]);
    kerns{4} = fft2(Gphi);

end