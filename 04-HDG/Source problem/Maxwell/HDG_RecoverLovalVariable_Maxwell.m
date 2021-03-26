function [wh,uh,ph] = HDG_RecoverLovalVariable_Maxwell(mymesh,...
        k,hat_variable,List_LocSol,List_LocSol_f)
    
    %% ------  Recover Local Solutions --------------------------
    Nw = (k+1)*(k+2)/2;
    Nu = 2*Nw;
    Np = Nw;
    N_local = Nw+Nu+Np;
    
    Nuhat_t = k+1;
    Nphat = k+1;
    
    N_global = Nuhat_t + Nphat;
    num_elements = mymesh.num_elements;
    %%% Local solvers
    %  ( qh )  =   Q*uhat + Qw*f
    %  ( uh )  =   U*uhat + Uw*f
    %
    % List_LocSol --> (Q,U)
    % List_LocSol_f-->(Qw*f,Uw*f)
    %%%%
    
    Result_matrix = zeros(N_local,num_elements, numeric_t);
    
    for ele_idx = 1:num_elements
        element_faces_list = mymesh.element_faces_list(ele_idx,:);
        % go through each face to do (Q,U)'*uhat
        
        for ll = 1:length(element_faces_list)
            face_idx = element_faces_list(ll);
            temp_hat = hat_variable ((face_idx-1)*N_global+1:face_idx*N_global);
            
            Result_matrix(:,ele_idx) = Result_matrix(:,ele_idx)...
                + List_LocSol(:,:,ele_idx,ll)*temp_hat;
            
        end
        % add (Qwf,Uwf)'
        Result_matrix(:,ele_idx) = Result_matrix(:,ele_idx)...
            +List_LocSol_f(:,ele_idx);
    end
    
    wh = Result_matrix(1:Nw,:);
    uh = Result_matrix(Nw+1:Nw+Nu,:);
    ph =  Result_matrix(Nw+Nu:end,:);
    
    
end