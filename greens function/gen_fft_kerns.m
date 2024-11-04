function kerns = gen_fft_kerns(kerns,sz,ind)
        
    GS = reshape(kerns{1},sz);
    Gphi = reshape(kerns{4},sz);

    hessGS = kerns{2};
    gradlapGS = kerns{3};

    hessGSxx = reshape(hessGS(:,:,1),sz);
    hessGSxy = reshape(hessGS(:,:,2),sz);
    hessGSyy = reshape(hessGS(:,:,3),sz);

    gradlapGSx = reshape(gradlapGS(:,:,1),sz);
    gradlapGSy = reshape(gradlapGS(:,:,2),sz);

    [zi,zj] = ind2sub(sz,ind);

    GS = circshift(GS,[zi,zj]);
    kerns{1} = fft2(GS);

    hessGSxx = circshift(hessGSxx,[zi,zj]);
    hessGSxx = fft2(hessGSxx);
    hessGSxy = circshift(hessGSxy,[zi,zj]);
    hessGSxy = fft2(hessGSxy);
    hessGSyy = circshift(hessGSyy,[zi,zj]);
    hessGSyy = fft2(hessGSyy);
    kerns{2} = cat(3,hessGSxx,hessGSxy,hessGSyy);

    gradlapGSx = circshift(gradlapGSx,[zi,zj]);
    gradlapGSx = fft2(gradlapGSx);
    gradlapGSy = circshift(gradlapGSy,[zi,zj]);
    gradlapGSy = fft2(gradlapGSy);
    kerns{3} = cat(3,gradlapGSx,gradlapGSy);

    Gphi = circshift(Gphi,[zi,zj]);
    kerns{4} = fft2(Gphi);


end