function err = L2Error_scalar_Square(uh_GQ_pts,...
        uexact, GQ_x, GQ_y, GQ1DRef_wts, hx, hy)
    
    [NGQ,num_ele] = size(GQ_x);
    err_list = zeros(num_ele,1,numeric_t);
    scale_factor = hx * hy * numeric_t('1/4.0');
    
    y_cut = 10;
    xpts = reshape(GQ_x,[],1);
    xpts = xpts';
    fpts = zeros(1, num_ele*NGQ, numeric_t);
    
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
       
        fpts(1,(ii-1)*NGQ+1:ii*NGQ) = uexact_pts(y_cut,:) - uh_GQ_pts(y_cut,:,ii);

    end
    err = sqrt(sum(err_list));
    err_list = sqrt(err_list);
    
    
    N_col = sqrt(num_ele);
    N_level = 1;
    stard_id = (N_level-1)*N_col*NGQ + 1;
    end_id = N_level*N_col*NGQ;
   
    
    figure
    plot(xpts(stard_id:end_id),fpts(stard_id:end_id),'*--')
    hold on
    yline(0,'r-');
    for j = 1:N_col-1
        xline(hx*j,'r--');
    end
    
    text = "y cut and GQ #" + num2str(y_cut) +" and level "+num2str(N_level);
    title(text)
    legend('u- uh at GQ')
    
    
end