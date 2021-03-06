function [uhat_out, uhatD_out] = ProjectUhat(mymesh, k_out, k_in, uhat, uD, uN, GQ1DRef_pts,GQ1DRef_wts)
    % project uhat from k_in to k_out
    % Since we use Legendre polynomial, the projection is simple injection
    % the only extra thing we need to handle is the boundary uD
    % Nueman bounary is not implemented
    
    Nuhat_in = k_in + 1;
    Nuhat_out = k_out + 1;
    num_elements = mymesh.num_elements;
    num_faces = mymesh.num_faces;
    
    uhat_out = zeros(Nuhat_out, num_faces, numeric_t);
    uhat_in = reshape(uhat, [Nuhat_in,num_faces]);
    
    uhat_out(1:Nuhat_in,:) = uhat_in(:,:);
    
    uhat_out = reshape(uhat_out,[Nuhat_out * num_faces, 1]);
    uhatD_out = numeric_t(sparse(Nuhat_out*num_faces,1));
    
    for element_idx = 1: num_elements
        
        ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);
        temp_element = mymesh.element_list(element_idx,:);
        vertice_list = mymesh.vertices_list(temp_element(:),:);
        
        uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
        
        [edge_len_list,~] = GetTriFaceInfo(vertice_list);
        
        for ii = 1:length(ele_face_idx_list)
            face_id = ele_face_idx_list(ii);
            bdry_flag = mymesh.f_type(face_id);
            
            if bdry_flag == 1 % dirichlet boundary 
                % dirichlet boundary face 
                % just simply enforce boundary conditions
                start_id=(face_id-1) * Nuhat_out+1;
                end_id = face_id * Nuhat_out;
                
                Jk = mymesh.Jacobian_list(element_idx);
                
                temp_proj = Project_F_to_Face(Jk,vertice_list,...
                    ii,uhat_dir_list(1,ii),edge_len_list(ii),...
                    k_out,uD,GQ1DRef_pts,GQ1DRef_wts);
                
                uhat_out(start_id:end_id,1) = temp_proj/(edge_len_list(ii)*numeric_t('0.5'));
                uhatD_out(start_id:end_id,1) = temp_proj/(edge_len_list(ii)*numeric_t('0.5'));
                
            end
        end
    end
    
end