% testing finite difference   

h = 0.1;
[X,Y] = meshgrid(-10:h:10);

% d2 - partial_{xx} (8th order)
d2 = zeros(9, 1);
d2(1) = -1/560;
d2(2) = 8/315;
d2(3) = -1/5;
d2(4) = 8/5;
d2(5) = -205/72;
d2(6) = 8/5;
d2(7) = -1/5;
d2(8) = 8/315;
d2(9) = -1/560;

% d3 - partial_{xxx} (6th order)
d3 = [-7/240; 3/10; -169/120; 61/30; 0; -61/30; 169/120; -3/10; 7/240];	

% d1 - partial_{x} (8th order)
d1 = [1/280; -4/105; 1/5; -4/5; 0; 4/5; -1/5; 4/105; -1/280];

gradlapx = zeros(9);
gradlapx(5,:) = d3.';
gradlapx = gradlapx + d2*d1.';
gradlapx = gradlapx / h^3;

q = X.^3.*Y.^2;

ind = find((X == 1) & (Y == 1));
[ii,jj] = ind2sub(size(X),ind);

subq = q(ii-4:ii+4,jj-4:jj+4);

val = sum(subq.*gradlapx,'all')

