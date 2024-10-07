% Finite difference test for gamma = 0
% zeroing gamma kills the nonlocal term
% making the equation easier to apply

h = 0.05;
xs = -10:h:10;
[X,Y] = meshgrid(xs);
gamma = 0;
beta = 5;
R = sqrt(X.^2 + Y.^2);
gval = green(R,beta,gamma);

% finite difference stencils
d1 = zeros(9,1);
d1(1) = 7/240;	
d1(2) = -2/5;
d1(3) = 169/60;
d1(4) = -122/15;
d1(5) = 91/8;
d1(6) = -122/15;
d1(7) = 169/60;
d1(8) = -2/5;
d1(9) = 7/240;
d1 = d1/h^4; 

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
d2 = d2 / h^2;

bilap = zeros(9,9);
bilap(:,5) = d1;
bilap(5,:) = bilap(5,:) + d1.';
bilap = bilap + 2*(d2*d2.');

sz = size(gval);
m = round(sz(1)/3);
n = round(2*sz(1)/3);
subf = gval(m-4:m+4,n-4:n+4);
err = abs(sum(bilap .* subf,'all') - beta*gval(m,n) ) / abs(gval(m,n))