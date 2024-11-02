%%%%%
%
% Continuous scattering with variable beta 
%
%%%%%

addpath(genpath('..'))

L = 50;
hs = 5./2.^(2:6); % (1:7)
errs = hs*0;

% Constant coefficients
a0 = 1;
b0 = 3; 
g0 = -1;

% Finding positive real roots
[rts,ejs] = find_roots(b0 / a0, g0 / a0);
k = rts((imag(rts) == 0) & (real(rts) > 0));

for ii = 1:numel(hs)

    % Parameters
    h = hs(ii);

    xs = -L:h:L;
    xl = -2*L:h:2*L;
    [~,n] = size(xs);
    [X,Y] = meshgrid(xs);
    [XL,YL] = meshgrid(xl);

    src = [0;0];
    targ = [XL(:).'; YL(:).'];

    % Perturbation to beta
    bbar = -0.5*X.*exp(-(X.^2 + Y.^2)/(2*(4*k)^2));
    beta = b0 + bbar;
    
    coefs = {a0,bbar};
    
    % RHS (Incident field)
    phiinc = exp(1i*k*X);
    rhs = k*bbar.*phiinc;
    rhs_vec = rhs(:);
    
    figure(1);
    tiledlayout(1,2);

    nexttile
    s = pcolor(X,Y,beta);
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
    [inds,corrs] = get_correct(rts,ejs,h);
    kerns = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),inds,corrs);

    ind = find((XL == 0) & (YL ==0));
    sz = size(XL);

    kerns = gen_fft_kerns(kerns,sz,ind);
    
    evalkerns = {kerns{1}, kerns{4}};
    
    % Solve with GMRES
    start = tic;
    mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,[],1e-12,200);
    mu = reshape(mu, size(X));
    t1 = toc(start);
    fprintf('%5.2e s : time to solve\n',t1)
    
    [phi, phi_n] = sol_eval_fft(mu,evalkerns);
    
    phi_tot = phiinc + phi;
    phi_n_tot = phi_n + k*phiinc;
    
   
    figure(2);
    tiledlayout(1,5)

    nexttile
    pc = pcolor(X,Y,real(mu));
    pc.EdgeColor = 'none';
    title('Re(\mu)')
    colorbar

    nexttile
    pc = pcolor(X,Y,real(phi_tot));
    pc.EdgeColor = 'none';
    title('Re(\phi)')
    colorbar

    nexttile
    pc = pcolor(X,Y,abs(phi_tot));
    pc.EdgeColor = 'none';
    title('|\phi|')
    colorbar

    nexttile
    pc = pcolor(X,Y,real(phi_n_tot));
    pc.EdgeColor = 'none';
    title('real(\phi_n)')
    colorbar

    nexttile
    pc = pcolor(X,Y,abs(phi_n_tot));
    pc.EdgeColor = 'none';
    title('|\phi_n|')
    colorbar
           
    % Calculate error with finite difference
    errs(ii) = get_fin_diff_err(X,Y,mu,phi_n_tot,phi_tot,X*0+a0,beta,g0,h);

end

%% Plotting errors 

figure(3)
loglog(hs,errs,'x-');

hold on
loglog(hs, 0.05*hs.^6, '--')

legend('error','h^6')

