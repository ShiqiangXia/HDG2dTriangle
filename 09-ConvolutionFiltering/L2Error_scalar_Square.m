function err = L2Error_scalar_Square(uh_GQ_pts,...
        uexact, GQ_x, GQ_y, GQ1DRef_wts, hx, hy)
    
    [NGQ,num_ele] = size(GQ_x);
    err_list = zeros(num_ele,1,numeric_t);
    scale_factor = hx * hy * numeric_t('1/4.0');
    
    for ii = 1:num_ele
        % get GQ points
        % total of NGQxNGQ points
        m = ones(NGQ,NGQ, numeric_t);
        x_list = m.*GQ_x(:,ii);
        x_list = reshape(x_list,[],1);
        
        y_list = m.*GQ_y(:,ii);
        y_list = reshape(y_list',[],1);
        
        uexact_pts = uexact([x_list,y_list]);
                
        uexact_pts = reshape(uexact_pts,[],NGQ);
        
        % make sure the matrix data matches with the structure of uh_GQ_pts
        % each row is the same y, x go from left to right
        uexact_pts = uexact_pts'; 
        
        
        diff2 = (uh_GQ_pts(:,:,ii) - uexact_pts).^2;
        
        % do quadratrue.
        err_list(ii,1) =...
           scale_factor*GQ1DRef_wts'*(diff2)*GQ1DRef_wts;

    end
    err = sqrt(sum(err_list));
    err_list = sqrt(err_list);
    
end