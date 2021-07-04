function [uh2star,qh2star] = HDG_2nd_Local_Postprocess_Elliptic(pb,mymesh,...
        GQ1DRef_pts,GQ1DRef_wts,...
        k_2star,k_star,tau,uhstar,source_f,uD,uN)
    
    % Idea: use uhstar as uhat in the local equation of higher order HDG
    % method and solve for uh* and qh*
    
    %% ---- step 1: Get local equations -----------------------------------
    
    [List_LocSol, List_LocSol_f, List_Ns]...
        = HDG_GetLocalEquations_Elliptic(pb, mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k_2star,tau,source_f);
    
    
    %% ------  Recover Local Solutions --------------------------
    
    Nu = (k_2star+1)*(k_2star+2)/2;
    Nq = 2*Nu;
    Nuhat = k_2star+1;
    dir_vec = GetDirVec(Nuhat); % correct the uhat oritation
    num_elements = mymesh.num_elements;
    
    V1D = Vandermonde1D(k_2star,GQ1DRef_pts);% legendre polynomail at GQ points
    
    Result_matrix = zeros(Nq+Nu,num_elements, numeric_t);
    
    for ele_idx = 1:num_elements
        
        element_faces_list = mymesh.element_faces_list(ele_idx,:);
        
        uhstar_coeff = uhstar(:,ele_idx);
        
        uhat_dir_list = mymesh.uhat_dir_list(ele_idx,:);

        % get uhat based on uhstar
        for ll = 1:length(element_faces_list)
            
            [face_r_list,face_s_list] = GetRefFaceQuadPts(ll,GQ1DRef_pts);
            [face_a_list,face_b_list] = RStoAB(face_r_list,face_s_list);
            
            V2Dstar_face = Vandermonde2D(k_star,face_a_list,face_b_list);
            
            uhstar_face_pts = V2Dstar_face*uhstar_coeff;
            
            temp_uhat = GQ1DRef_wts'*(uhstar_face_pts.*V1D);
            temp_uhat = temp_uhat';
            
            if uhat_dir_list(ll) == 0
                temp_uhat = temp_uhat.*dir_vec;
            end
            
            
            Result_matrix(:,ele_idx) = Result_matrix(:,ele_idx)...
                + List_LocSol(:,:,ele_idx,ll)*temp_uhat;
            
        end
        
        % add (Qwf,Uwf)'
        Result_matrix(:,ele_idx) = Result_matrix(:,ele_idx)...
            +List_LocSol_f(:,ele_idx);
        
        
        
    end
    
    qh2star = Result_matrix(1:Nq,:);
    uh2star = Result_matrix(Nq+1:end,:);
    
    
    
end