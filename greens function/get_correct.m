function [inds, corrs] = get_correct(rts,ejs,h,a0)
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

    tic

    inds = cell(1,4);

    A5 = [1 1 1 1 1;
        0 1 -1 0 0;
        0 0 0 1 -1;
        0 1 1 0 0;
        0 0 0 1 1];

    % kernmat(zi,zj) = kernmat(zi,zj) + tau(1);
    % kernmat(zi,zj+1) = kernmat(zi,zj+1) + tau(2);
    % kernmat(zi,zj-1) = kernmat(zi,zj-1) + tau(3);
    % kernmat(zi+1,zj) = kernmat(zi+1,zj) + tau(4);
    % kernmat(zi-1,zj) = kernmat(zi-1,zj) + tau(5);
    % kernmat(zi,zj+2) = kernmat(zi,zj+2) + tau(6);
    % kernmat(zi,zj-2) = kernmat(zi,zj-2) + tau(7);
    % kernmat(zi+2,zj) = kernmat(zi+2,zj) + tau(8);
    % kernmat(zi-2,zj) = kernmat(zi-2,zj) + tau(9);
    % kernmat(zi+1,zj+1) = kernmat(zi+1,zj+1) + tau(10);
    % kernmat(zi+1,zj-1) = kernmat(zi+1,zj-1) + tau(11);
    % kernmat(zi-1,zj+1) = kernmat(zi-1,zj+1) + tau(12);    
    % kernmat(zi-1,zj-1) = kernmat(zi-1,zj-1) + tau(13);

    A13 = [ones(1,13); ...
        0 1 -1 0 0 2 -2 0 0 1 -1 1 -1; ...
        0 0 0 1 -1 0 0 2 -2 1 1 -1 -1; ...
        0 1 1 0 0 4 4 0 0 1 1 1 1; ...
        0 0 0 1 1 0 0 4 4 1 1 1 1; ...
        zeros(1,9) 1 -1 -1 1; ...
        zeros(1,9) 1 1 -1 -1; ...
        zeros(1,9) 1 -1 1 -1; ...
        0 1 -1 0 0 8 -8 0 0 1 -1 1 -1; ...
        0 0 0 1 -1 0 0 8 -8 1 1 -1 -1; ...
        0 1 1 0 0 16 16 0 0 1 1 1 1; ...
        zeros(1,9) ones(1,4); ...
        0 0 0 1 1 0 0 16 16 1 1 1 1];   

    i1 = [0 0; 1 0; -1 0; 0 1; 0 -1];
    i2 = [0 0; 1 0; -1 0; 0 1; 0 -1; 2 0; -2 0; 0 2; 0 -2; 1 1; -1 1; 1 -1; -1 -1];
    
    inds{1} = i1;
    inds{2} = i2;
    inds{3} = i2;
    inds{4} = i1;
    
    % 1/2 r^2 log r^2 correction
    
    [~,z1] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    z1 = imag(z1)*1e12;
    [~,~,z2] = epstein_zeta(-4+1i*10^-12,1,0,1,1,0,0) ;
    z2 = imag(z2)*1e12;
    [~,~,z3] = epstein_zeta(-4+1i*10^-12,1,0,1,0,1,0) ;
    z3 = imag(z3)*1e12;

    b = [2*z1; 0; 0; z2/2+z3/8; z2/2+z3/8];
    b = b*h^4;
    tau = A5 \ b;
    c0 = 1/(4*pi*a0); % change this back to 1/(4 pi)

    valcor = c0*tau;

    % log(|r|^2) + 2 x^2 / r^2 + 1 
    
    [z0] = epstein_zeta(0+1i*10^-12,1,0,1) ;
    z0 = imag(z0)*1e12 ;
    [~,z1] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    z1 = imag(z1)*1e12;

    b = [z0 + log(h); 0; 0; z1; z1];
    b = b*h^2;
    tau0 = A5 \ b;

    [~,z1] = epstein_zeta(0+1i*10^-12,1,0,1,1,0,0) ;
    z1 = imag(z1)*1e12;
    [~,~,z2] = epstein_zeta(-2+1i*10^-12,1,0,1,0,1,0) ;
    z2 = imag(z2)*1e12;
    [~,~,z3] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    z3 = imag(z3)*1e12;

    b = [2*z1; 0; 0; 2*z3; 1/2*z2];
    b = b*h^2;
    tau1 = A5 \ b;

    hessxxcor = [c0*(2*tau0 + 2*tau1 ); zeros(8,1)];

    % 2*log(|r|) + 2 y^2 / r^2 + 1 

    b = [2*z1; 0; 0; 1/2*z2; 2*z3];
    b = b*h^2;
    tau1 = A5 \ b;

    hessyycor = [c0*(2*tau0 + 2*tau1 ); zeros(8,1)];

    % 2*x*y / r^2

    [~,~,z2] = epstein_zeta(-2+1i*10^-12,1,0,1,0,1,0) ;
    z2 = imag(z2)*1e12;

    b = [0; 0; 0; 0; 0; 1/2*z2; 0; 0; 0; 0; 0; 0; 0];
    b = b*h^2;
    tau = A13 \ b;

    hessxycor = 2*c0*tau;

    % 4 x / r^2

    [~,z1] = epstein_zeta(0+1i*10^-12,1,0,1,1,0,0) ;
    z1 = imag(z1)*1e12;
    [~,~,z2] = epstein_zeta(-2+1i*10^-12,1,0,1,0,1,0) ;
    z2 = imag(z2)*1e12;
    [~,~,z3] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    z3 = imag(z3)*1e12;

    b = [0; 2*z1; 0; 0; 0; 0; 0; 1/2*z2; 2*z3; 0; 0; 0; 0];
    b = b*h;
    tau = A13 \ b;

    gradlapxcor = 4*c0*tau;

    % 4 y / r^2

    b = [0; 0; 2*z1; 0; 0; 0; 1/2*z2; 0; 0; 2*z3; 0; 0; 0];
    b = b*h;
    tau = A13 \ b;

    gradlapycor = 4*c0*tau;

    A13 = [ones(1,13); ...
        0 1 -1 0 0 2 -2 0 0 1 -1 1 -1; ...
        0 0 0 1 -1 0 0 2 -2 1 1 -1 -1; ...
        0 1 1 0 0 4 4 0 0 1 1 1 1; ...
        0 0 0 1 1 0 0 4 4 1 1 1 1; ...
        zeros(1,9) 1 -1 -1 1; ...
        zeros(1,9) 1 1 -1 -1; ...
        zeros(1,9) 1 -1 1 -1; ...
        0 1 -1 0 0 8 -8 0 0 1 -1 1 -1; ...
        0 0 0 1 -1 0 0 8 -8 1 1 -1 -1; ...
        0 1 1 0 0 16 16 0 0 1 1 1 1; ...
        zeros(1,9) ones(1,4); ...
        0 0 0 1 1 0 0 16 16 1 1 1 1]; 

    % |r|^3 correction

    [z0,~] = epstein_zeta_int(-3,1,0,1,0,0,0) ;
    [~,z1] = epstein_zeta_int(-5,1,0,1,1,0,0) ;

    b = [-z0; 0; 0; -2/5*z1; -2/5*z1];
    b = b*h^5;
    tau = A5 \ b;
    c0 = 1/4/4/gamma(1+3/2)^2; 

    phicor = c0*tau;

    toc



    hesscor = cat(3,hessxxcor,hessxycor,hessyycor);
    gradlapcor = cat(3,gradlapxcor,gradlapycor);

    corrs = {valcor, hesscor, gradlapcor, phicor};



end