function Explicit_Functional_Error_Terms(func_type,pde_ype,para,...
        mymesh,uh,qh,uhat,vh,ph,vhat,GQ1DRef_pts,GQ1DRef_wts,k,tau,post_flag)
    
    if strcmp(func_type,'1')
        Nu = (k+1)*(k+2)/2;
        Nuhat  = k+1;
        NGQ = length(GQ1DRef_pts);
        num_elements = mymesh.num_elements;
        dir_vec = GetDirVec(Nuhat); % correct the uhat oritation 
        dir_vec = dir_vec';
        % Get Gauss Quadpoints on the square
        [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts);
        % Map Gauss Quadpoints to the reference triangle
        [r_list,s_list] = ABtoRS(a_list,b_list);
        
        if strcmp(pde_ype,'1') % Poission
            % need qexact, pexact
            [qexact_1,qexact_2,pexact_1,pexact_2]= MyParaParse(para,...
                'qexact_1','qexact_2','pexact_1','pexact_2');
            if post_flag == 0
                V2D = Vandermonde2D(k,a_list,b_list);
                [V2Dr,V2Ds] = GradVandermonde2D(k,a_list,b_list);
            elseif post_flag == 1
                error('Post-processed has not implemented yet');
            end
            V1D = Vandermonde1D(k,GQ1DRef_pts);% legendre polynomail at GQ points
            % compute the error terms over all elements
            Err1_elewise_list = zeros(num_elements,1,numeric_t);
            Err2_elewise_list = zeros(num_elements,1,numeric_t);
            Err3_elewise_list = zeros(num_elements,1,numeric_t);
            Err3_elewise_list = zeros(num_elements,1,numeric_t);
            
            for element_idx = 1: num_elements
                
                temp_element = mymesh.element_list(element_idx,:);
                vertice_list = mymesh.vertices_list(temp_element(:),:);
                Jk = mymesh.Jacobian_list(element_idx);
                [~,Inv_AffineMap] = Ref_Tri_Map(Jk,vertice_list);
                % Gauss points on any element
                [x_list,y_list] = RStoXY(r_list,s_list,Jk,vertice_list);
                
                uh_coeff = uh(:,element_idx);
                vh_coeff = vh(:,element_idx);
                qh_coeff_1 = qh(1:Nu,element_idx);
                qh_coeff_2 = qh(Nu+1:end,element_idx);
                ph_coeff_1 = ph(1:Nu,element_idx);
                ph_coeff_2 = ph(Nu+1:end,element_idx);
                
                
                % qexact, pexact at those Gauss points
                q1_VD = qexact_1([x_list,y_list]);q1_VD = reshape(q1_VD,[],NGQ);
                
                q2_VD = qexact_2([x_list,y_list]);q2_VD = reshape(q2_VD,[],NGQ);
                
                p1_VD = pexact_1([x_list,y_list]);p1_VD = reshape(p1_VD,[],NGQ);
                
                p2_VD = pexact_2([x_list,y_list]);p2_VD = reshape(p2_VD,[],NGQ);
                
                % qh,ph at Gauss points
                
                qh_pts1 = V2D * (qh_coeff_1);qh_pts1 = reshape(qh_pts1,[],NGQ);
                
                qh_pts2 = V2D * (qh_coeff_2);qh_pts2 = reshape(qh_pts2,[],NGQ);
                
                ph_pts1 = V2D * (ph_coeff_1);ph_pts1 = reshape(ph_pts1,[],NGQ);
                
                ph_pts2 = V2D * (ph_coeff_2);ph_pts2 = reshape(ph_pts2,[],NGQ);
         
            
            % err_1 (q-qh,p-ph)
                temp_formula_1 = ((q1_VD-qh_pts1).*(p1_VD-ph_pts1) ...
                                +  (q2_VD-qh_pts2).*(p2_VD-ph_pts2) );
                            
                Err1_elewise_list(element_idx,1) =...
                    Jk*GQ1DRef_wts'*(temp_formula_1.*Jacobian_rs_to_ab )*GQ1DRef_wts;
            
            
            %d/dx = d/dr * dr/dx + d/ds * ds/dx
            grad_vh_1 = V2Dr * vh_coeff * Inv_AffineMap(1,1)...
                        + V2Ds * vh_coeff * Inv_AffineMap(2,1); grad_vh_1 = reshape(grad_vh_1,[],NGQ);
            
            grad_vh_2 = V2Dr * vh_coeff * Inv_AffineMap(1,2)...
                        + V2Ds * vh_coeff * Inv_AffineMap(2,2); grad_vh_2 = reshape(grad_vh_2,[],NGQ);
            
            grad_uh_1 = V2Dr * uh_coeff * Inv_AffineMap(1,1)...
                        + V2Ds * uh_coeff * Inv_AffineMap(2,1); grad_uh_1 = reshape(grad_uh_1,[],NGQ);
            
            grad_uh_2 = V2Dr * uh_coeff * Inv_AffineMap(1,2)...
                        + V2Ds * uh_coeff * Inv_AffineMap(2,2); grad_uh_2 = reshape(grad_uh_2,[],NGQ);
            
            % err_2 (q-qh,ph+grad(vh)) + (qh+grad(uh),p-ph)
            
            temp_formula_2 = ...
                 (q1_VD-qh_pts1).*(ph_pts1 +grad_vh_1)...
                + (q2_VD-qh_pts2).*(ph_pts2 +grad_vh_2)...
                + (qh_pts1 + grad_uh_1).*(p1_VD-ph_pts1)...
                + (qh_pts2 + grad_uh_2).*(p2_VD-ph_pts2);
            
            Err2_elewise_list(element_idx,1) = ...
                Jk*GQ1DRef_wts'*(temp_formula_2.*Jacobian_rs_to_ab )*GQ1DRef_wts;
            
            % err_3 <qhat*n - q*n, vh-vhat> + < uh- uhat, phat*n-p*n>
            ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);
            uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
            [e_list,n1,n2,n3] = GetTriFaceInfo(vertice_list);
            
            for ii = 1:length(ele_face_idx_list)
                face_id = ele_face_idx_list(ii);
                start_id=(face_id-1)*Nuhat+1;
                end_id = face_id*Nuhat;
                
                uhat_coeff = uhat(start_id:end_id);
                vhat_coeff = vhat(start_id:end_id);
                
                
                
            end
            
            
            
            
            end
            
            
            
        end
    else
        error ('This type of error terms calculation has not been implemented yet.')
    end
    
        
end