function [qh_Neig,uh_Neig] = HDG_Eig_RecoverLovalVariable(mymesh,...
        k,lam_list,uhat_Neig,List_LocSol,List_LocSol_f,Neig)
    
    Nu = (k+1)*(k+2)/2;
    Nq = 2*Nu;
    Nuhat = k+1;
    num_elements = mymesh.num_elements;
    
    Result_matrix = zeros(Nq+Nu,num_elements, Neig,numeric_t);
    
    List_Qw = List_LocSol_f(1:Nq,:,:);
    List_Uw = List_LocSol_f(Nq+1:end,:,:);
    
    List_Q = List_LocSol(1:Nq,:,:,:);
    List_U = List_LocSol(Nq+1:end,:,:,:);
    

    
    for ii=1:Neig
        
        lam = lam_list(1,ii);
        uh_hat = uhat_Neig(:,ii);
        
        for ele_idx = 1:num_elements
            element_faces_list = mymesh.element_faces_list(ele_idx,:);
            
            Uw = List_Uw(:,:,ele_idx);
            
            temp_mat = eye(Nu,Nu,numeric_t) - lam * Uw;
            
            for ll = 1:length(element_faces_list)
                face_idx = element_faces_list(ll);
                temp_uhat = uh_hat ((face_idx-1)*Nuhat+1:face_idx*Nuhat);
                
                Result_matrix(Nq+1:end,ele_idx,ii) = Result_matrix(Nq+1:end,ele_idx,ii)...
                + (temp_mat\List_U(:,:,ele_idx,ll))*temp_uhat;
            
                Result_matrix(1:Nq,ele_idx,ii) = Result_matrix(1:Nq,ele_idx,ii)...
                + (List_Q(:,:,ele_idx,ll))*temp_uhat;
            end
            
            Result_matrix(1:Nq,ele_idx,ii) = Result_matrix(1:Nq,ele_idx,ii)...
                + lam*List_Qw(:,:,ele_idx)*Result_matrix(Nq+1:end,ele_idx,ii);
            
        end
    end
    
    qh_Neig = Result_matrix(1:Nq,:,:);
    
    uh_Neig = Result_matrix(Nq+1:end,:,:);
    
end