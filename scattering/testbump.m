
h = 0.1;
[X,Y] = meshgrid(-10:h:10);
xl = -2;
xr = 6;
yl = -3;
yr = 8;

[psi,psix,psiy,psixx,psixy,psiyy] = bumpfunc(X,Y,xl,xr,yl,yr);

figure(1);
s = pcolor(X,Y,psi);
s.EdgeColor = 'none';
colorbar


figure(2);
s = pcolor(X,Y,psiy);
s.EdgeColor = 'none';
colorbar


figure(3);
s = pcolor(X,Y,psixy);
s.EdgeColor = 'none';
colorbar

