function err1 = get_fin_diff_err_helm(X,Y,utot,h,coefs,xloc,yloc,zk)

    [~,ind] = min((X(:) - xloc(:)).^2 + (Y(:) - yloc(:)).^2);
    [ii, jj] = ind2sub(size(X),ind);
    % disp(phi(ii,jj))
    
    % d2 - partial_{xx} (8th order)
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

    lap = zeros(9);
    lap(5,:) = d2.';
    lap(:,5) = lap(:,5) + d2;
    lap = lap / h^2;

    V = coefs{1};

    % Residual error of total solution 
    usub = utot(ii-4:ii+4,jj-4:jj+4);

    err1 = abs(sum(lap.*usub,'all') + zk^2*(1 + V(ii,jj))*utot(ii,jj)) ;
    
end