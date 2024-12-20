% checking the derivatives for jeremy's ridges

load('cracks.mat');
fout = fout;
foutx = foutx;
fouty = fouty;
foutxx = foutxx;
foutxy = foutxy;
foutyy = foutyy;


ind = find(fout == max(fout(:)));
[i,j] = ind2sub(size(fout),ind);

% stencils

d1 = [1/280; -4/105; 1/5; -4/5; 0; 4/5; -1/5; 4/105; -1/280];

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

(fouty(i,j) - fout(i-4:i+4,j).'*d1) / fouty(i,j)
(foutx(i,j) - fout(i,j-4:j+4)*d1) / foutx(i,j)

(foutxx(i,j) - fout(i,j-4:j+4)*d2) / foutxx(i,j)
(foutyy(i,j) - fout(i-4:i+4,j).'*d2) / foutyy(i,j)


dxy = d1*d1.';

(foutxy(i,j) - sum(dxy.*fout(i-4:i+4,j-4:j+4),'all')) / foutxy(i,j)