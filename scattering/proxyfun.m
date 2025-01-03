% TODO: change the documentation as multipple chunkers are adopted... 

function [Kpxy,nbr] = proxyfun(slf,nbr,l,ctr,x,h,rts,ejs,V,pr,pw,t1,t2)
%PROXYFUN proxy function utility for kernels defined on chunkers
%
% We assume that the kernel for a given source and target is an 
% opdims(1) x opdims(2) matrix. Matrix entry indices use the 
% following convention: Similarly index j corresponds 
% to the (j-1)/opdims(2) + 1 boundary point for slf points and 
% index i corresponds to the (i-1)/opdims(1) + 1 boundary point
% for the nbr points (from the point of view of this routine slf 
% points will be sources where nbr points are potential targets)
%
% this function primarily calls the provided kernel routine
% but makes some attempts at efficiency (kernel calls corresponding to
% the same point are not repeated).
%
% ~ tree and index related inputs ~
% slf - set of relevant indices for self box (not necessarily equivalent
%       to the appropriate indices within the chunker object)
% nbr - set of indices in neighbor boxes (to check for inclusion inside 
%        proxy surface)
% l - box side length at current level
% ctr - center of current box
%
% ~ other inputs ~
% chnkrs - chunker object or an array of chunker objects
% kern - should be a kernel function of the form 
%             submat = kern(src, targ, srctau, targtau)
%        or multiple kernel functions for the multiple chunker case
% opdims_mat - dimensions of operator, opdims(1) output dim, opdims(2) input
%            dim of operator (distinct from spatial dim)
%            size = (2,nchunker,nchunker)
% pr - a (dim,_) array of proxy points (dim is spatial dim here)
% ptau - a (dim,_) array of proxy unit tangents
% pw - is a set of smooth integration weights for the proxy points
% pin - function handle, takes in rescaled and recentered points in the 
%       nbr array and determines if they're within the proxy surface
% l2scale - boolean type that determines if we should 
%    rescale the matrix by l2scale. the default value is false. 
% 

% scaled proxy points and weights (no scaling necessary on tangents)

lmax = max(l);
pxy = pr*lmax + ctr(:);
pw = lmax^2*pw;
nbr = nbr(sum(((x(:,nbr) - ctr)/lmax).^2) < 1.5^2);

%%%%%

src  = pxy;
targ = x(:,slf);
vals = green(src,targ,rts,ejs);

v1 = vals{1};
v2a= squeeze(vals{2}(:,:,1));
v2b= squeeze(vals{2}(:,:,2));
v2c= squeeze(vals{2}(:,:,3));
v3a= squeeze(vals{3}(:,:,1));
v3b= squeeze(vals{3}(:,:,2));

v4 = vals{4};

v2lap = v2a + v2c;

%%%%%%%

a0 = V{1};
abar = V{2}(slf);
b0 = V{3};
bbar = V{4}(slf);
g0 = V{5};
gbar = V{6};
alphax = V{7}(slf);
alphay = V{8}(slf);
alphaxx = V{9}(slf);
alphaxy = V{10}(slf);
alphayy = V{11}(slf);
nu = V{end};
alphalap = alphaxx + alphayy;

xmat = diag(alphax)*v3a + diag(alphay)*v3b + ...
    0.5*diag(alphalap)*v2lap + ...
    0.5*(1-nu)*(2*diag(alphaxy)*v2b-diag(alphayy)*v2a-diag(alphaxx)*v2c)-...
    0.5*diag((a0*bbar-abar*b0)./a0)*v1 - diag(0.5*abar.*g0./a0)*v4;
xmat = xmat*diag(pw);


vals = green(targ,src,rts,ejs);
ymat = vals{1}.'*h^2;
zmat = vals{4}.'*h^2;
%sum(abs(xmat(:)))/sum(abs(ymat(:)))
%sum(abs(xmat(:)))/sum(abs(zmat(:)))
%ymat = ymat*sum(abs(xmat(:)))/sum(abs(ymat(:)));
%zmat = zmat*sum(abs(xmat(:)))/sum(abs(zmat(:)));
ymat = t1*ymat;
zmat = t2*zmat;
Kpxy = [xmat.';ymat.';zmat.'];
%Kpxy = xmat.';
%%%%%%%%%%%%%%%%%%


end
