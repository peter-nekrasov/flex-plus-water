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
% format of cell arrays mirrors that of green.m


    inds = cell(1,4);
    corrs = cell(1,4);
    
    inds{1} = [0 0; 0 1; 0 -1; 1 0; -1 0];
    inds{4} = [0 0; 0 1; 0 -1; 1 0; -1 0];
    
    % r^2 log r corrections
    
    [~,z1] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    z1 = imag(z1)*1e12;
    [~,~,z2] = epstein_zeta(-4+1i*10^-12,1,0,1,1,0,0) ;
    z2 = imag(z2)*1e12;
    [~,~,z3] = epstein_zeta(-4+1i*10^-12,1,0,1,0,1,0) ;
    z3 = imag(z3)*1e12;
    
    A = [1 1 1 1 1;
        0 1 -1 0 0;
        0 0 0 1 -1;
        0 1 1 0 0;
        0 0 0 1 1];
    b = [2*z1; 0; 0; z2/2+z3/8; z2/2+z3/8];
    b = b*h^4;
    tau = A \ b;
    c0 = 1/(8*pi);
    corrs{1} = c0*tau;
    corrs{4} = c0*tau;

end