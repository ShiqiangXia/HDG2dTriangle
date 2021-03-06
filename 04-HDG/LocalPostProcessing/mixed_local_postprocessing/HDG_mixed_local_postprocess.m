function uh_out = HDG_mixed_local_postprocess(mymesh,uh_ave, qh_RT, k_in, k_out, GQ1DRef_pts,GQ1DRef_wts)
    % the idea is as follows
    % use qh* (in RT_k) from HDG local processing to postprocessing u variable
    % (grad uh**, grad w) = -(qh*, grad w), where w in P_{k+1}
    %  (uh**, 1)   = (uh, 1)
    
    num_elements = mymesh.num_elements;
    NRT = (k_in + 1) * (k_in + 3); % dimension of RT_kin
    Nuh_out = (k_out + 2) * (k_out +1) /2; % dimension of P_{k_out}
    uh_out = zeros(Nuh_out, num_elements, numeric_t);
    
    % step 1: obtain local matrix on ref element
    % (ur,ur), (ur, us),(us,us) for P__{k+2}
    % (RTk_1, ur), (RTk_1, us),(RTk_2, ur), (RTk_2, us)
    
    [RT1_ur, RT1_us, RT2_ur, RT2_us] = ...
        Volume_Int_RT_du(k_in, k_out, GQ1DRef_pts,GQ1DRef_wts);
    
    [Aurur,Aurus,Ausus] = Volume_Int_du_du(k_out,GQ1DRef_pts,GQ1DRef_wts);
    
    % step2: Go through each element 
    % build local matrix and do postprocessing
    
    for element_idx = 1: num_elements
        temp_element = mymesh.element_list(element_idx,:);
        
        %element_faces_list  = mymesh.element_faces_list(element_idx,:);
        
        vertice_list = mymesh.vertices_list(temp_element(:),:);
        Jk = mymesh.Jacobian_list(element_idx);
        
        %uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
        
        qh_RT_coef = qh_RT(:,element_idx);
        
        W_mat = PostScalar_Assemble(Jk,vertice_list,Aurur,Aurus,Ausus);
        W_RHS = MixedPostScalar_GetRHS(Jk,vertice_list,RT1_ur, RT1_us, RT2_ur, RT2_us,qh_RT_coef, uh_ave(1,element_idx));
        
        uh_out(:,element_idx) = W_mat\W_RHS;
        
    end
    
end