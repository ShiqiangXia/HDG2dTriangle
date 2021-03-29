function [M0,N0] = HDG_Eig_AssembleMlambda_Elliptic...
        (lam,mymesh,k,List_Ns,List_LocSol,List_LocSol_f)
    
    Nu = (k+1)*(k+2)/2;
    Nq = 2*Nu;
    Nuhat = k+1;
    
    List_Uw = List_LocSol_f(Nq+1:end,:,:); % Nu x Nu
    List_U = List_LocSol(Nq+1:end,:,:,:);
    
    num_faces = mymesh.num_faces;
    num_elements = mymesh.num_elements;
    
    M0 = numeric_t(sparse(Nuhat*num_faces,Nuhat*num_faces));
    N0 = numeric_t(sparse(Nuhat*num_faces,Nuhat*num_faces));
    
    for element_idx = 1: num_elements
        ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);

        Uw = List_Uw(:,:,element_idx);
        temp_mat = eye(Nu,Nu,numeric_t) - lam * Uw;
        
        for ii = 1:length(ele_face_idx_list)
            face_id = ele_face_idx_list(ii);
          
            start_id=(face_id-1)*Nuhat+1;
            end_id = face_id*Nuhat;
         
            bdry_flag = mymesh.f_type(face_id);
            
            if bdry_flag == 0
                
                Bd_Int_mat = List_Ns(:,:,element_idx,ii)';
                S1 = - Bd_Int_mat * List_LocSol_f(:,:,element_idx);
                
                for jj = 1:length(ele_face_idx_list)
                    temp_id = ele_face_idx_list(jj);
                    temp_start = (temp_id-1)*Nuhat+1;
                    temp_end = temp_id*Nuhat;
                    
                    Loc_U = List_U(:,:,element_idx,jj);
                    
                    
                    S2 = (temp_mat\Loc_U);
                    S3 = (temp_mat\Uw);
                    
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