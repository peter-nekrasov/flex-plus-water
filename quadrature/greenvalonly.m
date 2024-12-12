function val = greenvalonly(targ,x,y,rts,ejs)

    src = [x(:).'; y(:).'];
    gf = green(src,targ,rts,ejs);
    val = reshape(gf{1},size(x));
    

end