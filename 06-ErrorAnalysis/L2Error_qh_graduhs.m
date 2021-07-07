function [err,err_list] = L2Error_qh_graduhs(mymesh, k1, qh1, k2, uh2, GQ1DRef_pts,GQ1DRef_wts)
    % compute the error of qh - (-grad uh) (they may have different degrees)
    
    Nu1 = (k1+2)*(k1+1)/2;
    
            
    qh1_1 = qh1(1:Nu1,:);
    qh1_2 = qh1(Nu1+1:end,:);
    
    NGQ = length(GQ1DRef_pts);
    num_elements = mymesh.num_elements;
    
    err_list = zeros(num_elements,1,numeric_t);
    [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts);
    
    V2D_1 = Vandermonde2D(k1,a_list,b_list);  % scalar Pk1
    [V2D_r, V2D_s] = GradVandermonde2D(k2,a_list,b_list);
    for element_idx = 1: num_elements
        temp_element = mymesh.element_list(element_idx,:);
        vertice_list = mymesh.vertices_list(temp_element(:),:);
        Jk = mymesh.Jacobian_list(element_idx);
        
        qh1_pts1 = V2D_1 * (qh1_1(:,element_idx));
        qh1_pts1 = reshape(qh1_pts1,[],NGQ);
        qh1_pts2 = V2D_1 * (qh1_2(:,element_idx));
        qh1_pts2 = reshape(qh1_pts2,[],NGQ);
        
        qh2_pts1 = - GetGradUhPts(1,Jk,vertice_list,V2D_r,V2D_s,uh2(:,element_idx));
        qh2_pts1 = reshape(qh2_pts1,[],NGQ);
        qh2_pts2 = - GetGradUhPts(1,Jk,vertice_list,V2D_r,V2D_s,uh2(:,element_idx));
        qh2_pts2 = reshape(qh2_pts2,[],NGQ);
        
        diff1 = (qh1_pts1 - qh2_pts1).^2 ;
        
        diff2 = (qh1_pts2 - qh2_pts2).^2 ;
        
        temp_err1 = Jk*GQ1DRef_wts'*(diff1.*Jacobian_rs_to_ab)*GQ1DRef_wts;
        temp_err2 = Jk*GQ1DRef_wts'*(diff2.*Jacobian_rs_to_ab)*GQ1DRef_wts;
        
        err_list(element_idx,1) = temp_err1 + temp_err2;
 
    end
    
    err = sqrt(sum(err_list));
    err_list = sqrt(err_list);
    
end