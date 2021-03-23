function [Global_b]= HDG_AssembleGlobalEquation_RHS(mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,uD,uN, List_LocSol_f,List_Ns)
    
    %% ----- Assemble Global Matrix and right hand side --------
    
    Nuhat = k+1;
    
    num_elements = mymesh.num_elements;
    num_faces = mymesh.num_faces;
    % Preparation
    Global_b = numeric_t(sparse(Nuhat*num_faces,1));
    
    
    for element_idx = 1: num_elements
        
        ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);
        temp_element = mymesh.element_list(element_idx,:);
        vertice_list = mymesh.vertices_list(temp_element(:),:);
        
        uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
        
        [edge_len_list,~] = GetTriFaceInfo(vertice_list);
        for ii = 1:length(ele_face_idx_list)
            face_id = ele_face_idx_list(ii);
          
            start_id=(face_id-1)*Nuhat+1;
            end_id = face_id*Nuhat;
         
            bdry_flag = mymesh.f_type(face_id);
            
            if bdry_flag == 0 % interior face
                
                Bd_Int_mat = List_Ns(:,:,element_idx,ii)';
  
                Global_b(start_id:end_id,1) =Global_b(start_id:end_id,1)-Bd_Int_mat*List_LocSol_f(:,element_idx);
                
                
                     
            elseif bdry_flag == 1 % dirichlet boundary 
                % dirichlet boundary face 
                % just simply enforce boundary conditions
                
                Jk = mymesh.Jacobian_list(element_idx);

                Global_b(start_id:end_id,1) = Project_F_to_Face(Jk,vertice_list,...
                    ii,uhat_dir_list(1,ii),edge_len_list(ii),...
                    k,uD,GQ1DRef_pts,GQ1DRef_wts);
                
                
            elseif bdry_flag == 2 % neuman boundary
                % Will implement later
                error(' Boundary type not implemented yet.')
            else 
                error('Wrong boundary type.')
            end
            
        end
    end
    
end