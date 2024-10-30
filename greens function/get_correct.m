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
    
    % c0 r^2 log r + c1 |r|^3 corrections
    
    [~,z1] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    z1 = imag(z1)*1e12;
    [~,~,z2] = epstein_zeta(-4+1i*10^-12,1,0,1,1,0,0) ;
    z2 = imag(z2)*1e12;
    [~,~,z3] = epstein_zeta(-4+1i*10^-12,1,0,1,0,1,0) ;
    z3 = imag(z3)*1e12;
    
    A5 = [1 1 1 1 1;
        0 1 -1 0 0;
        0 0 0 1 -1;
        0 1 1 0 0;
        0 0 0 1 1];

    b = [2*z1; 0; 0; z2/2+z3/8; z2/2+z3/8];
    b = b*h^4;
    tau = A5 \ b;
    c0 = 1/(8*pi);

    corrs{1} = c0*tau;


    [z0,~] = epstein_zeta_int(-3,1,0,1,0,0,0) ;
    [~,z1] = epstein_zeta_int(-5,1,0,1,1,0,0) ;

    b = [-z0; 0; 0; -2/5*z1; -2/5*z1];
    b = b*h^5;
    tau = A5 \ b;
    c0 = 1/4/8/gamma(1+3/2)^2; 

    corrs{4} = c0*tau;

end