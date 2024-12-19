%%%%%
%
% Testing convergence of the adjointed Lippman-Schwinger equation 
% using finite differences
%
%%%%%

clear 
close all
addpath(genpath('..'))

L = 6;
hs = (1/2).^(0:8);
hsmall = hs(end-1);
hsmallest = hs(end);
errsinf = hs*0;
errs2 = hs*0;
errconv1 = hs*0;
errconv2 = hs*0;
errconvinf = hs*0;

muconv = load('mu.mat').mu;

[Xerr,Yerr] = meshgrid(-1:0.25:1);


for ii = 1:numel(hs)

    h = hs(ii);
    
    xs = -L/2:h:L/2;
    xl = -L:h:L;
    [~,n] = size(xs);
    [X,Y] = meshgrid(xs);
    [XL,YL] = meshgrid(xl);
    [N,~] = size(X);
    
    [coefs, H] = bump3(X,Y,3,0.5,100); % remove gbar from coefs vector
    
    a0 = coefs{1}; 
    b0 = coefs{3}; 
    g0 = coefs{5}; 
    
    % Finding positive real roots
    [rts,ejs] = find_roots(b0 / a0, g0 / a0);
    k = rts((imag(rts) == 0) & (real(rts) > 0));
    ejs = ejs/a0;
    
    src = [0;0];
    targ = [XL(:).'; YL(:).'];
    
    % RHS (Incident field)
    k1 = 2*sqrt(2)*k/3;
    k2 = k/3;
    phiinc = exp(1i*k1*X+1i*k2*Y);
    [rhs_vec, rhs] = get_rhs_vec(coefs,k1,k2,phiinc);
    
    figure(1);
    tiledlayout(1,3);
    
    nexttile
    s = pcolor(X,Y,coefs{1} + coefs{2});
    s.EdgeColor = 'None';
    colorbar
    title('\alpha')
    drawnow
    
    nexttile
    s = pcolor(X,Y,coefs{2} + coefs{3});
    s.EdgeColor = 'None';
    colorbar
    title('\beta')
    drawnow
    
    nexttile
    s = pcolor(X,Y,real(rhs));
    s.EdgeColor = 'None';
    colorbar
    title('rhs')
    drawnow
    
    
    % Constructing integral operators
    [inds,corrs] = get_correct(h,a0);
    kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs), h, inds,corrs);
    
    ind = find((XL == 0) & (YL ==0));
    sz = size(XL);
    
    kerns = gen_fft_kerns(kerns,sz,ind);

    evalkerns = {kerns{1}, kerns{4}};
        
    % Solve with GMRES
    start = tic;
    mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,10,1e-15,200);
    mu = reshape(mu, size(X));
    t1 = toc(start);
    fprintf('%5.2e s : time to solve\n',t1)

    [phi_interp, phi_n_interp] = sol_eval_fft(mu,evalkerns);

    phi_tot = phi_interp + phiinc;
    phi_n_tot = phi_n_interp + k*phiinc;


    if h > hsmall     % Interpolating the density

        muhat = fftshift(fft2(mu));
        px = -L/2:hsmallest:L/2+h-hsmallest;
        [Xsmallest,Ysmallest] = meshgrid(px);
        Np = length(px);
        npad = (Np-N)/2;
        muhat = padarray2(muhat,npad);
        muhat = ifftshift(muhat);
        mu_interp = ifft2(muhat)*Np^2/N^2;

        tmpinterp = mu_interp(1:1537,1:1537);
        errconv2(ii) = sqrt(sum(abs(tmpinterp(:) - muconv(:)).^2)) / sqrt(sum(abs(muconv(:)).^2));
        errconv1(ii) = (sum(abs(tmpinterp(:) - muconv(:)))) / (sum(abs(muconv(:))));
        errconvinf(ii) = max(abs(tmpinterp(:) - muconv(:))) / max(abs(muconv(:)));

        muhat = fftshift(fft2(mu));
        px = -L/2:hsmall:L/2+h-hsmall;
        [Xsmall,Ysmall] = meshgrid(px);
        Np = length(px);
        npad = (Np-N)/2;
        muhat = padarray2(muhat,npad);
        muhat = ifftshift(muhat);
        mu_interp = ifft2(muhat)*Np^2/N^2;

        % Evaluating solution using interpolated density
    
        pxl = -L-h+hsmall:hsmall:L+h-hsmall;
        [XLsmall,YLsmall] = meshgrid(pxl);
        
        src = [0;0];
        targ = [XLsmall(:).'; YLsmall(:).'];
    
        [inds,corrs] = get_correct(hsmall,a0);
        kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs), hsmall, inds,corrs);
    
        ind = find((XLsmall == 0) & (YLsmall ==0));
        sz = size(XLsmall);
    
        kerns = gen_fft_kerns(kerns,sz,ind);
    
        evalkerns = {kerns{1}, kerns{4}};
    
        [phi_interp, phi_n_interp] = sol_eval_fft(mu_interp,evalkerns);
    
        phiinc = exp(1i*k1*Xsmall+1i*k2*Ysmall);
        phi_tot_interp = phi_interp + phiinc;
        phi_n_tot_interp = phi_n_interp + k*phiinc;

    else

        muhat = fftshift(fft2(mu));
        px = -L/2:hsmallest:L/2+h-hsmallest;
        [Xsmallest,Ysmallest] = meshgrid(px);
        Np = length(px);
        npad = (Np-N)/2;
        muhat = padarray2(muhat,npad);
        muhat = ifftshift(muhat);
        mu_interp = ifft2(muhat)*Np^2/N^2;

        tmpinterp = mu_interp(1:1537,1:1537);
        errconv2(ii) = sqrt(sum(abs(tmpinterp(:) - muconv(:)).^2)) / sqrt(sum(abs(muconv(:)).^2));
        errconv1(ii) = (sum(abs(tmpinterp(:) - muconv(:)))) / (sum(abs(muconv(:))));
        errconvinf(ii) = max(abs(tmpinterp(:) - muconv(:))) / max(abs(muconv(:)));
    
        px = -L/2:hsmall:L/2+h-hsmall;
        [Xsmall,Ysmall] = meshgrid(px);

        phi_n_tot_interp = phi_n_tot(ismember(X,Xsmall) & ismember(Y,Ysmall));
        phi_n_tot_interp = reshape(phi_n_tot_interp, size(Xsmall));
        phi_tot_interp = phi_tot(ismember(X,Xsmall) & ismember(Y,Ysmall));
        phi_tot_interp = reshape(phi_tot_interp, size(Xsmall));
        mu_interp = mu(ismember(X,Xsmall) & ismember(Y,Ysmall));
        mu_interp = reshape(mu_interp, size(Xsmall));

    end



    figure(2);
    tiledlayout(2,2)
    
    nexttile
    pc = pcolor(X,Y,real(mu));
    pc.EdgeColor = 'none';
    title('Re(\mu) original')
    colorbar

    nexttile
    pc = pcolor(X,Y,real(phi_tot));
    pc.EdgeColor = 'none';
    title('Re(\phi) original')
    colorbar
    
    nexttile
    pc = pcolor(Xsmall,Ysmall,real(mu_interp));
    pc.EdgeColor = 'none';
    title('Re(\mu) interpolated')
    colorbar
    
    nexttile
    pc = pcolor(Xsmall,Ysmall,real(phi_tot_interp));
    pc.EdgeColor = 'none';
    title('Re(\phi) interpolated')
    colorbar
    
    % figure(3);
    % tiledlayout(1,4)
    % 
    % 
    % nexttile
    % pc = pcolor(Xsmall,Ysmall,real(phi_tot));
    % pc.EdgeColor = 'none';
    % title('Re(\phi)')
    % colorbar
    % 
    % nexttile
    % pc = pcolor(Xsmall,Ysmall,abs(phi_tot));
    % pc.EdgeColor = 'none';
    % title('|\phi|')
    % colorbar
    % 
    % nexttile
    % pc = pcolor(Xsmall,Ysmall,real(phi_n_tot));
    % pc.EdgeColor = 'none';
    % title('real(\phi_n)')
    % colorbar
    % 
    % nexttile
    % pc = pcolor(Xsmall,Ysmall,abs(phi_n_tot));
    % pc.EdgeColor = 'none';
    % title('|\phi_n|')
    % colorbar
           

    [coefs, H] = bump3(Xsmall,Ysmall,3,0.5,100);
    

    % Calculate error with finite difference
    err = get_fin_diff_err3(Xsmall,Ysmall,mu_interp,phi_n_tot_interp,phi_tot_interp,hsmall,coefs,Xerr(:),Yerr(:));
    
    errsinf(ii) = max(err)
    errs2(ii) = mean(err)

end

%% Plotting
 
figure(4);


nexttile
loglog(hs,errs2,'o-','DisplayName','$\ell_2$ - finite difference')
hold on
loglog(hs,errsinf,'o-','DisplayName','$\ell_\infty$ - finite difference')
hold on
loglog(hs,errconv2,'o-','DisplayName','$L_2$ - self convergence')
hold on

loglog(hs(3:end-1),2000*hs(3:end-1).^6,'k--','DisplayName','$h^6$')
hold on


xlim([0.5*min(hs), 2*max(hs)])
%ylim([0.1*min(errs2(1:end-1)) 10*max(errs2(1:end-1))])
xlabel('h')
ylabel('Relative error')
title('Residual')
legend('Interpreter','latex','Location','eastoutside')
