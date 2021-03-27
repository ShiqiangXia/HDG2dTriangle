function [Global_b]= HDG_AssembleGlobalEquation_RHS_Maxwell(mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,uhat_t_D, List_LocSol_f,List_Ns)
    
    %% ----- Assemble Global Matrix and right hand side --------
    
    Nuhat_t = k+1;
    Nphat =  k+1;
    N_global =Nuhat_t + Nphat;
    
    num_elements = mymesh.num_elements;
    num_faces = mymesh.num_faces;
    % Preparation
    Global_b = numeric_t(sparse(N_global*num_faces,1));
    
    
    for element_idx = 1: num_elements
        
        ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);
        temp_element = mymesh.element_list(element_idx,:);
        vertice_list = mymesh.vertices_list(temp_element(:),:);
        
        uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
        
        [edge_len_list,~] = GetTriFaceInfo(vertice_list);
        for ii = 1:length(ele_face_idx_list)
            face_id = ele_face_idx_list(ii);
          
            start_id=(face_id-1)*N_global+1;
            end_id = face_id*N_global;
         
            bdry_flag = mymesh.f_type(face_id);
            
            if bdry_flag == 0 % interior face
                
                Bd_Int_mat = List_Ns(:,:,element_idx,ii)';
  
                Global_b(start_id:end_id,1) =Global_b(start_id:end_id,1)-Bd_Int_mat*List_LocSol_f(:,element_idx);
                
                
                     
            elseif bdry_flag == 1 % dirichlet boundary 
                % dirichlet boundary face 
                % just simply enforce boundary conditions
                
                Jk = mymesh.Jacobian_list(element_idx);
                
                % set uhat_t = uhat_t_D on border 
                %     phat = 0
                Global_b(start_id:start_id+Nuhat_t-1,1) = Project_F_to_Face(Jk,vertice_list,...
                    ii,uhat_dir_list(1,ii),edge_len_list(ii),...
                    k,uhat_t_D,GQ1DRef_pts,GQ1DRef_wts);
                
                
            elseif bdry_flag == 2 % neuman boundary
                % Will implement later
                error(' Boundary type not implemented yet.')
            else 
                error('Wrong boundary type.')
            end
            
        end
    end
    
end