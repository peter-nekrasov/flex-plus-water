function val = phivalonly(targ,x,y,rts,ejs)

    src = [x(:).'; y(:).'];
    gf = green(src,targ,rts,ejs);
    val = reshape(gf{4},size(x));
    

end