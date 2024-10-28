function [inds, corrs] = get_correct(rts,ejs,h)
% Getting corrections for the kernels in the Lippman-Schwinger eq
%
% input:
% - rts: float vector - roots of polynomial
% - ejs: float vector - coefficients in partial fraction expansion
% - h: float - grid spacing
%
% output: 
% - inds: cell array - indices of corrections
% - corrs: cell array - corresponding corrections
%
% format of cell array: 
% inds - 

    [~,z1] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    z1 = imag(z1)*1e12;
    [~,~,z2] = epstein_zeta(-4+1i*10^-12,1,0,1,1,0,0) ;
    z2 = imag(z2)*1e12;
    [~,~,z3] = epstein_zeta(-4+1i*10^-12,1,0,1,0,1,0) ;
    z3 = imag(z3)*1e12;

    [rts, ejs] = find_roots(beta,gamma);

    %%% PERFORM DIAGONAL CORRECTIONS HERE
    % NEED TO KNOW CONSTANT PART 

    [zi,zj] = ind2sub(size(kern_struct{1}),ind);

    A = [1 1 1 1 1;
        0 1 -1 0 0;
        0 0 0 1 -1;
        0 1 1 0 0;
        0 0 0 1 1];
    b = [2*z1; 0; 0; z2/2+z3/8; z2/2+z3/8];
    b = b*h^4;
    tau = A \ b;
    c0 = sum(ejs.*rts.^4)/(8*pi);

    kernmat(zi,zj) = kernmat(zi,zj) + c0*tau(1);
    kernmat(zi,zj+1) = kernmat(zi,zj+1) + c0*tau(2);
    kernmat(zi,zj-1) = kernmat(zi,zj-1) + c0*tau(3);
    kernmat(zi+1,zj) = kernmat(zi+1,zj) + c0*tau(4);
    kernmat(zi-1,zj) = kernmat(zi-1,zj) + c0*tau(5);



end