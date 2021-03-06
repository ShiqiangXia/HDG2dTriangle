function [Global_M]= HDG_AssembleGlobalEquation_LHS_Maxwell(mymesh,...
        k,tau_t,tau_n,epsilon, List_LocSol,List_Ns)
    
    %% ----- Assemble Global Matrix and right hand side --------
    
    Nuhat_t = k+1;
    Nphat =  k+1;
    N_global =Nuhat_t + Nphat;
    
    num_elements = mymesh.num_elements;
    num_faces = mymesh.num_faces;
    % Preparation
    Global_M = numeric_t(sparse(N_global*num_faces,N_global*num_faces));
    
    Id_mtrix = eye(Nuhat_t,Nuhat_t, numeric_t);
    
    
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
                
                % put matrix block at the right position
                %<qh*n+tau*uh,mu>
                
                for jj = 1:length(ele_face_idx_list)
                    temp_id = ele_face_idx_list(jj);
                    temp_start = (temp_id-1)*N_global+1;
                    temp_end = temp_id*N_global;
                    
                    Global_M(start_id:end_id,temp_start:temp_end) =...
                        Global_M(start_id:end_id,temp_start:temp_end)+...
                       Bd_Int_mat * List_LocSol(:,:,element_idx,jj) ;
                end

                
                
                temp_Id_mat = blkdiag(tau_t*Id_mtrix,tau_n*epsilon*Id_mtrix);
                
                
                Global_M(start_id:end_id,start_id:end_id) = ...
                    Global_M(start_id:end_id,start_id:end_id) - temp_Id_mat*edge_len_list(ii)*0.5;
                     
            elseif bdry_flag == 1 % dirichlet boundary 
                % dirichlet boundary face 
                % just simply enforce boundary conditions
                temp_Id_mat = blkdiag(Id_mtrix,Id_mtrix);
                Global_M(start_id:end_id,start_id:end_id)= temp_Id_mat*edge_len_list(ii)*0.5;
                
          
            elseif bdry_flag == 2 % neuman boundary
                % Will implement later
                error(' Boundary type not implemented yet.')
            else 
                error('Wrong boundary type.')
            end
            
        end
    end
    
end