function val = hessxyonly(targ,x,y,rts,ejs)

    src = [x(:).'; y(:).'];
    gf = green(src,targ,rts,ejs);
    gf = gf{2};
    val = gf(:,:,2);
    val = reshape(val,size(x));

end