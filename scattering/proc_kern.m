function kern_struct = proc_kern(kern_struct,h,ind,beta,gamma)

    val = kern_struct{1};
    hessxx = kern_struct{2};
    hessxy = kern_struct{3};
    hessyy = kern_struct{4};
    gradlapx = kern_struct{5};
    gradlapy = kern_struct{6};
    phi = kern_struct{7};

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

    for ii = 1:numel(kern_struct)
        kernmat = kern_struct{ii};
        if ii == 1
                A = [1 1 1 1 1;
                    0 1 -1 0 0;
                    0 0 0 1 -1;
                    0 1 1 0 0;
                    0 0 0 1 1];
                b = [2*z1; 0; 0; z2/2+z3/8; z2/2+z3/8];
                b = b*h^4;
                tau = A \ b;
                c0 = 1; %sum(ejs.*rts.^4)/(8*pi);
            
                kernmat(zi,zj) = kernmat(zi,zj) + c0*tau(1);
                kernmat(zi,zj+1) = kernmat(zi,zj+1) + c0*tau(2);
                kernmat(zi,zj-1) = kernmat(zi,zj-1) + c0*tau(3);
                kernmat(zi+1,zj) = kernmat(zi+1,zj) + c0*tau(4);
                kernmat(zi-1,zj) = kernmat(zi-1,zj) + c0*tau(5);
        end
        kernhat = circshift(kernmat,[zi,zj]);
        kern_struct{ii} = fft2(kernhat)*h*h;
    end


end