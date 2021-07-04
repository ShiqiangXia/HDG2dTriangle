function [uhat_out, uhatD_out] = ...
        ProjectUh_LocalPostprocessed(mymesh, k_out, k_in,...
        uhstar, uD, uN, GQ1DRef_pts,GQ1DRef_wts)
    
    
    Nuhat_out = k_out + 1;
    num_elements = mymesh.num_elements;
    num_faces = mymesh.num_faces;
    COEFF  = numeric_t('0.5');
    uhat_out = zeros(Nuhat_out*num_faces, 1, numeric_t);
    
    dir_vec = GetDirVec(Nuhat_out); % correct the uhat oritation
    V1D = Vandermonde1D(k_out,GQ1DRef_pts);% legendre polynomail at GQ points
    
    for element_idx = 1: num_elements
        ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);
        temp_element = mymesh.element_list(element_idx,:);
        vertice_list = mymesh.vertices_list(temp_element(:),:);
        uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
        [edge_len_list,~] = GetTriFaceInfo(vertice_list);
        uhstar_coeff = uhstar(:,element_idx);
        for ii = 1:length(ele_face_idx_list)
            face_id = ele_face_idx_list(ii);
            start_id=(face_id-1) * Nuhat_out+1;
            end_id = face_id * Nuhat_out;
            bdry_flag = mymesh.f_type(face_id);
            if bdry_flag == 1 % dirichlet boundary 
                % dirichlet boundary face 
                % just simply enforce boundary conditions
                
                Jk = mymesh.Jacobian_list(element_idx);
                temp_proj = Project_F_to_Face(Jk,vertice_list,...
                    ii,uhat_dir_list(1,ii),edge_len_list(ii),...
                    k_out,uD,GQ1DRef_pts,GQ1DRef_wts);
                
                uhat_out(start_id:end_id,1) = temp_proj/(edge_len_list(ii)*0.5);
                uhatD_out(start_id:end_id,1) = temp_proj/(edge_len_list(ii)*0.5);
                
            elseif bdry_flag == 0 % interior face
                
                [face_r_list,face_s_list] = GetRefFaceQuadPts(ii,GQ1DRef_pts);
                [face_a_list,face_b_list] = RStoAB(face_r_list,face_s_list);
            
                V2Dstar_face = Vandermonde2D(k_in,face_a_list,face_b_list);
            
                uhstar_face_pts = V2Dstar_face*uhstar_coeff;
                
                temp_uhat = GQ1DRef_wts'*(uhstar_face_pts.*V1D);
                temp_uhat = temp_uhat';
            
                if uhat_dir_list(ii) == 0
                    temp_uhat = temp_uhat.*dir_vec;
                end
                uhat_out(start_id:end_id,1) = uhat_out(start_id:end_id,1)...
                    + COEFF  * temp_uhat ;
                
            
                
            end
        end
    end
    
end
    