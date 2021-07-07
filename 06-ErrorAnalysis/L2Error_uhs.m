function [err,err_list]=L2Error_uhs(mymesh, k1, uh1, k2, uh2, GQ1DRef_pts,GQ1DRef_wts,flag)
    % compute error of  uh1 - uh2 (they may have different degrees)
    % compute the L2 error of ||uh1-uh2||
    % get basic info
    NGQ = length(GQ1DRef_pts);
    num_elements = mymesh.num_elements;
    
    err_list = zeros(num_elements,1,numeric_t);
    [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts);
    
    if flag == 0
        V2D_1 = Vandermonde2D(k1,a_list,b_list);  % scalar Pk1
    elseif flag == 1
        [V2D_1,~] = RTVandermonde2D(k1,a_list,b_list); % 1st comp. of RT_k
    elseif flag == 2
        [~,V2D_1] = RTVandermonde2D(k1,a_list,b_list);% 2nd comp. of RT_k
    end
    
    V2D_2 = Vandermonde2D(k2,a_list,b_list);  % scalar Pk2
    
    for element_idx = 1: num_elements
        Jk = mymesh.Jacobian_list(element_idx);
        
        uh_pts1 = V2D_1 * (uh1(:,element_idx));
        uh_pts1 = reshape(uh_pts1,[],NGQ);
        uh_pts2 = V2D_2 * (uh2(:,element_idx));
        uh_pts2 = reshape(uh_pts2,[],NGQ);
        
        diff2 = (uh_pts1 - uh_pts2).^2;
        
        % do quadratrue.
        err_list(element_idx,1) =...
            Jk*GQ1DRef_wts'*(diff2.*Jacobian_rs_to_ab)*GQ1DRef_wts;
    end
    err = sqrt(sum(err_list));
    err_list = sqrt(err_list);

end   