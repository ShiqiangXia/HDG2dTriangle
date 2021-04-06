function HDG_Local_Postprocess_Elliptic(mymesh,k,uh,qh,uhat,tau,GQ1DRef_pts,GQ1DRef_wts)
    
    num_elements = mymesh.num_elements;
    Nqstar = (k+1)*(k+3);   % RT_k
    Nustar = (k+2)*(k+3)/2; % P_{k+1}
    
    
    [A_qtau,Buuhat3,RTuuhat3] = RT_LocalMatrix(k,GQ1DRef_pts,GQ1DRef_wts);
    
    
    qstar = zeros(Nqstar,num_elements,numeric_t);
    ustar = zeros(Nustar, num_elements, numeric_t);
    
    % go through each element and do local post-processing
    
    for element_idx = 1: num_elements
        temp_element = mymesh.element_list(element_idx,:);
        
        element_faces_list  = mymesh.element_faces_list(element_idx,:);
        
        vertice_list = mymesh.vertices_list(temp_element(:),:);
        Jk = mymesh.Jacobian_list(element_idx);
        
        uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
        
        uh_coef = uh(:,element_idx);
        qh_coef = qh(:,element_idx);
        
        RT_mat = RT_Assemble(k,Jk,vertice_list, A_qtau,Buuhat3,RTuuhat3);
        
        RT_RHS = RT_GetRHS(k,Jk,vertice_list,uhat_dir_list,A_qtau,Buuhat3,uh_coef,qh_coef,uhat,element_faces_list,tau);
        
        qstar(:,element_idx) = RT_mat\RT_RHS;
        
        W_mat = ;
        W_RHS = ;
        
        ustar(:,element_idx) = W_mat\W_RHS;
        
    end
    
end
    