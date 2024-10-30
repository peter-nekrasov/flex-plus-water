function err = get_fin_diff_err(X,Y,mu,phi_n_tot,phi_tot,alpha,beta,gamma,h)


    ind = intersect(find(X == 5), find(Y == 5));
    [ii, jj] = ind2sub(size(X),ind);
    % disp(phi(ii,jj))
    
    % finite difference stencils 
    % d1 - partial_{xxxx} (6th order) % improve order?
    d1 = zeros(9,1);
    d1(1) = 7/240;	
    d1(2) = -2/5;
    d1(3) = 169/60;
    d1(4) = -122/15;
    d1(5) = 91/8;
    d1(6) = -122/15;
    d1(7) = 169/60;
    d1(8) = -2/5;
    d1(9) = 7/240;
    
    % d2 - partial_{xxyy} (8th order)
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
    
    bilap = zeros(9,9);
    bilap(:,5) = d1;
    bilap(5,:) = bilap(5,:) + d1.';
    bilap = bilap + 2*(d2*d2.');
    bilap = bilap / h^4;
    
    % Error of scattered part
    % phi_n_sub = phi_n(ii-4:ii+4,jj-4:jj+4);
    % first = a0*sum(bilap.*phi_n_sub,'all') ;
    % second = -beta(ii,jj)*phi_n(ii,jj);
    % third = g0*phi(ii,jj);
    % rhs = k*bbar(ii,jj)*phiinc(ii,jj);
    % err = abs(first + second + third - rhs) / max(abs(phi_n_tot(:)))
    
    % Residual error of total solution 
    phi_n_tot_sub = phi_n_tot(ii-4:ii+4,jj-4:jj+4);
    first = alpha*sum(bilap.*phi_n_tot_sub,'all') ;
    second = -beta(ii,jj).*phi_n_tot(ii,jj);
    third = gamma*phi_tot(ii,jj);
    err = abs(first + second + third); % / abs(sum(h^2*mu(:)));

end