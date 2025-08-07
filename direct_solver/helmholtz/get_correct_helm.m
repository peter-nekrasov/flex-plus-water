function [inds, corr] = get_correct_helm(h)
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

    inds = cell(1);
    corr = cell(1);

    A5 = [1 1 1 1 1;
        0 1 -1 0 0;
        0 0 0 1 -1;
        0 1 1 0 0;
        0 0 0 1 1];

    i1 = [0 0; 1 0; -1 0; 0 1; 0 -1];
    
    inds{1} = i1;
    
    c0 = 1/(4*pi); % change this back to 1/(4 pi)

    % log(|r|^2) 
    
    [z0] = epstein_zeta(0+1i*10^-12,1,0,1) ;
    z0 = imag(z0)*1e12 ;
    [~,z1] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    z1 = imag(z1)*1e12;

    b = [z0 + log(h); 0; 0; z1; z1];
    tau0 = A5 \ b;

    corr{1} = c0*(2*tau0 );

end