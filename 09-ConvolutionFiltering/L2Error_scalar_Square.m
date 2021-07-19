function err = L2Error_scalar_Square(uh_GQ_pts,...
        uexact_GQ_pts, GQ1DRef_wts, hx, hy)
    
    [~,~,num_ele] = size(uh_GQ_pts);
    err_list = zeros(num_ele,1,numeric_t);
    scale_factor = hx * hy * numeric_t('1/4.0');
    
    
    for ii = 1:num_ele
        
        diff2 = (uh_GQ_pts(:,:,ii) - uexact_GQ_pts(:,:,ii)).^2;
        
        % do quadratrue.
        err_list(ii,1) =...
           scale_factor*GQ1DRef_wts'*(diff2)*GQ1DRef_wts;
       
    end
    err = sqrt(sum(err_list));
    err_list = sqrt(err_list);
    
end