function [wh_Neig,uh_Neig,ph_Neig] = HDG_Eig_RecoverLovalVariable_Maxwell(mymesh,...
        k,lam_list,hat_var_Neig,List_LocSol,List_LocSol_f,Neig)
    
    Nw = (k+1)*(k+2)/2;
    Nu = 2*Nw;
    Np = Nw;
    N_local = Nw+Nu+Np;
    
    Nuhat_t = k+1;
    Nphat = k+1;
    
    N_global = Nuhat_t + Nphat;
    
    num_elements = mymesh.num_elements;
    
    Result_matrix = zeros(N_local,num_elements, Neig,numeric_t);
    
    List_Wj = List_LocSol_f(1:Nw,:,:);
    List_Uj = List_LocSol_f(Nw+1:Nw+Nu,:,:);
    List_Pj = List_LocSol_f(Nw+Nu+1:end,:,:);
    
    List_W = List_LocSol(1:Nw,:,:,:);
    List_U = List_LocSol(Nw+1:Nw+Nu,:,:,:);
    List_P = List_LocSol(Nw+Nu+1:end,:,:,:);
    

    
    for ii=1:Neig
        
        lam = lam_list(1,ii);
        hat_var = hat_var_Neig(:,ii);
        
        for ele_idx = 1:num_elements
            element_faces_list = mymesh.element_faces_list(ele_idx,:);
            
            Uj = List_Uj(:,:,ele_idx);
            
            temp_mat = eye(Nu,Nu,numeric_t) - lam * Uj;
            
            for ll = 1:length(element_faces_list)
                face_idx = element_faces_list(ll);
                temp_hat = hat_var ((face_idx-1)*N_global+1:face_idx*N_global);
                
                % (Id - lam*Uj) uh = U*hat
                % uh = (Id-lam*Uh)^-1 * U * uhat
                Result_matrix(Nw+1:Nw+Nu,ele_idx,ii) = Result_matrix(Nw+1:Nw+Nu,ele_idx,ii)...
                + (temp_mat\List_U(:,:,ele_idx,ll))*temp_hat;
                % wh = W*hat + Wj*lam*uh
                % 1. W*hat 
                Result_matrix(1:Nw,ele_idx,ii) = Result_matrix(1:Nw,ele_idx,ii)...
                + (List_W(:,:,ele_idx,ll))*temp_hat;
                % ph = P*hat + Pj*lam*uh
                % 1. P*hat
                Result_matrix(Nw+Nu+1:end,ele_idx,ii) = Result_matrix(Nw+Nu+1:end,ele_idx,ii)...
                + (List_P(:,:,ele_idx,ll))*temp_hat;
            end
            % 2. Wj*lam*uh
            Result_matrix(1:Nw,ele_idx,ii) = Result_matrix(1:Nw,ele_idx,ii)...
                + lam*List_Wj(:,:,ele_idx)*Result_matrix(Nw+1:Nw+Nu,ele_idx,ii);
            % 2. Pj*lam*uh
            Result_matrix(Nw+Nu+1:end,ele_idx,ii) = Result_matrix(Nw+Nu+1:end,ele_idx,ii)...
                + lam*List_Pj(:,:,ele_idx)*Result_matrix(Nw+1:Nw+Nu,ele_idx,ii);
            
        end
    end
    
    wh_Neig = Result_matrix(1:Nw,:,:);
    
    uh_Neig = Result_matrix(Nw+1:Nw+Nu,:,:);
    
    ph_Neig = Result_matrix(Nw+Nu:end,:,:);
    
    
end