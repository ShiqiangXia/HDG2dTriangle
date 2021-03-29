function [M0,N0] = HDG_Eig_AssembleMlambda_Maxwell...
        (lam,mymesh,k,List_Ns,List_LocSol,List_LocSol_f)
    
    Nw = (k+1)*(k+2)/2;
    Nu = 2*Nw;
    Np = Nw;
    N_local = Nw+Nu+Np;
    
    Nuhat_t = k+1;
    Nphat = k+1;
    
    N_global = Nuhat_t + Nphat;
    
    List_Uj = List_LocSol_f(Nw+1:Nw+Nu,:,:); % Nu x Nu
    
    List_U = List_LocSol(Nw+1:Nw+Nu,:,:,:);
    
    num_faces = mymesh.num_faces;
    num_elements = mymesh.num_elements;
    
    M0 = numeric_t(sparse(N_global*num_faces,N_global*num_faces));
    N0 = numeric_t(sparse(N_global*num_faces,N_global*num_faces));
    
    for element_idx = 1: num_elements
        ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);

        Uj = List_Uj(:,:,element_idx);
        temp_mat = eye(Nu,Nu,numeric_t) - lam * Uj;
        
        for ii = 1:length(ele_face_idx_list)
            face_id = ele_face_idx_list(ii);
          
            start_id=(face_id-1)*N_global+1;
            end_id = face_id*N_global;
         
            bdry_flag = mymesh.f_type(face_id);
            
            if bdry_flag == 0
                
                Bd_Int_mat = List_Ns(:,:,element_idx,ii)';
                S1 = - Bd_Int_mat * List_LocSol_f(:,:,element_idx);
                
                for jj = 1:length(ele_face_idx_list)
                    temp_id = ele_face_idx_list(jj);
                    temp_start = (temp_id-1)*N_global+1;
                    temp_end = temp_id*N_global;
                    
                    Loc_U = List_U(:,:,element_idx,jj);
                    
                    
                    S2 = (temp_mat\Loc_U);
                    S3 = (temp_mat\Uj);
                    
                    M0(start_id:end_id,temp_start:temp_end) =...
                        M0(start_id:end_id,temp_start:temp_end)...
                       +S1*(S2);
                   
                    N0(start_id:end_id,temp_start:temp_end) =...
                        N0(start_id:end_id,temp_start:temp_end)...
                       +S1*S3*S2;
                   
                end
                 
            elseif bdry_flag == 1 % dirichlet boundary 
                % do nothing
                continue;
            elseif bdry_flag == 2 % neuman boundary
                % Will implement later
                error(' Boundary type not implemented yet.')
            else 
                error('Wrong boundary type.')
            end
        end
        
    end
    
    N0 = lam*N0 + M0;
    
end