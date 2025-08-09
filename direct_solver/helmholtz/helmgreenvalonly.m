function val = helmgreenvalonly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = helm2d.green_cell_helm(zk,src,targ);
    val = reshape(gf{1},size(x));
    

end