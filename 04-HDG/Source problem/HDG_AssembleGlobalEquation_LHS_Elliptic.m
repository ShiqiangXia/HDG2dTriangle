function [Global_M]= HDG_AssembleGlobalEquation_LHS_Elliptic(mymesh,...
        k,tau, List_LocSol,List_Ns)
    
    %% ----- Assemble Global Matrix and right hand side --------
    
    Nuhat = k+1;
    
    num_elements = mymesh.num_elements;
    num_faces = mymesh.num_faces;
    % Preparation
    Global_M = numeric_t(sparse(Nuhat*num_faces,Nuhat*num_faces));
    Id_mtrix = eye(Nuhat,Nuhat, numeric_t);
    
    for element_idx = 1: num_elements
        
        ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);
        temp_element = mymesh.element_list(element_idx,:);
        vertice_list = mymesh.vertices_list(temp_element(:),:);

        [edge_len_list,~] = GetTriFaceInfo(vertice_list);
        for ii = 1:length(ele_face_idx_list)
            face_id = ele_face_idx_list(ii);
          
            start_id=(face_id-1)*Nuhat+1;
            end_id = face_id*Nuhat;
         
            bdry_flag = mymesh.f_type(face_id);
            
            if bdry_flag == 0 % interior face
                
                Bd_Int_mat = List_Ns(:,:,element_idx,ii)';
                
                % put matrix block at the right position
                %<qh*n+tau*uh,mu>
                for jj = 1:length(ele_face_idx_list)
                    temp_id = ele_face_idx_list(jj);
                    temp_start = (temp_id-1)*Nuhat+1;
                    temp_end = temp_id*Nuhat;
                    
                    Global_M(start_id:end_id,temp_start:temp_end) =...
                        Global_M(start_id:end_id,temp_start:temp_end)+...
                       Bd_Int_mat * List_LocSol(:,:,element_idx,jj) ;
                end

                % -<tau uhat,mu>
                Global_M(start_id:end_id,start_id:end_id) = ...
                    Global_M(start_id:end_id,start_id:end_id) - tau*Id_mtrix*edge_len_list(ii)*0.5;
                     
            elseif bdry_flag == 1 % dirichlet boundary 
                % dirichlet boundary face 
                % just simply enforce boundary conditions
              
                Global_M(start_id:end_id,start_id:end_id)= Id_mtrix*edge_len_list(ii)*0.5;
                
          
            elseif bdry_flag == 2 % neuman boundary
                % Will implement later
                error(' Boundary type not implemented yet.')
            else 
                error('Wrong boundary type.')
            end
            
        end
    end
    
end