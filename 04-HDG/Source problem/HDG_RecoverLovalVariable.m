function [qh,uh] = HDG_RecoverLovalVariable(mymesh,...
        k,uh_hat,List_LocSol,List_LocSol_f)
    
    %% ------  Recover Local Solutions --------------------------
    Nu = (k+1)*(k+2)/2;
    Nq = 2*Nu;
    Nuhat = k+1;
    num_elements = mymesh.num_elements;
    %%% Local solvers
    %  ( qh )  =   Q*uhat + Qw*f
    %  ( uh )  =   U*uhat + Uw*f
    %
    % List_LocSol --> (Q,U)
    % List_LocSol_f-->(Qw*f,Uw*f)
    %%%%
    
    Result_matrix = zeros(Nq+Nu,num_elements, numeric_t);
    
    for ele_idx = 1:num_elements
        element_faces_list = mymesh.element_faces_list(ele_idx,:);
        % go through each face to do (Q,U)'*uhat
        for ll = 1:length(element_faces_list)
            face_idx = element_faces_list(ll);
            temp_uhat = uh_hat ((face_idx-1)*Nuhat+1:face_idx*Nuhat);
            
            Result_matrix(:,ele_idx) = Result_matrix(:,ele_idx)...
                + List_LocSol(:,:,ele_idx,ll)*temp_uhat;
            
        end
        % add (Qwf,Uwf)'
        Result_matrix(:,ele_idx) = Result_matrix(:,ele_idx)...
            +List_LocSol_f(:,ele_idx);
    end
    
    qh = Result_matrix(1:Nq,:);
    uh = Result_matrix(Nq+1:end,:);
    
    
end