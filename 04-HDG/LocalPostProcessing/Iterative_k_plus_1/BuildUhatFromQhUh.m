function uhat = BuildUhatFromQhUh(mymesh, k,qh,uh,tau, uhatD, List_Ns)
    
    % go through each element and build uhat from qh and uh
    % use the transmission condition
    %        formula uhat = tau+/T uh+  + tau-/T + uh-  + 1/T [[qh]]
    %         here T = tau+ + tau-,  [[qh]] = qh+ * n+  + qh- * n-
    
    % to get the right coefficients: multiply both side by a test func mu
    % we get 
    % <uhat, mu> 
    % = ( <qh+ * n + tau+ * uh+, mu> + <qh- * n + tau- * uh-, mu> )/(tau+ + tau-)
    
    
    Nuhat = k + 1 ;
    num_faces = mymesh.num_faces;
    COEFF =  numeric_t('1.0')/(tau+ tau);
    
    uhat = zeros(Nuhat * num_faces, 1, numeric_t);
    
    for element_idx = 1: num_elements
        
        ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);
        temp_element = mymesh.element_list(element_idx,:);
        vertice_list = mymesh.vertices_list(temp_element(:),:);
        [edge_len_list,~] = GetTriFaceInfo(vertice_list);
        
        temp_vec = [qh(:,element_idx); uh(:,element_idx)];
        
        for ii = 1:length(ele_face_idx_list)
            face_id = ele_face_idx_list(ii);
          
            start_id=(face_id-1)*Nuhat+1;
            end_id = face_id*Nuhat;
         
            bdry_flag = mymesh.f_type(face_id);
            
            scale_factor = numeric_t('1.0')/(edge_len_list(ii)*0.5);
            
            if bdry_flag == 0 % interior face
                
                Bd_Int_mat = List_Ns(:,:,element_idx,ii)'; % transpose
                
                uhat(start_id:end_id,1) = uhat(start_id:end_id,1)...
                    + COEFF * Bd_Int_mat * temp_vec * scale_factor ;
                
            elseif bdry_flag == 1 % dirichlet boundary 
                uhat(start_id:end_id,1) = uhatD(start_id:end_id,1) * scale_factor;
            end
        end
        
    end
    
    
    
    
end