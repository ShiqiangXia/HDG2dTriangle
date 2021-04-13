function [Est_elewise_list,Est1_elewise_list,...
        Est2_elewise_list,Est3_elewise_list,...
        Est4_elewise_list,Est5_elewise_list] ...
                = Functional_Eh_Estimate_Terms_By_Extraction(func_type,pde_ype,mymesh,...
                uhstar,qhstar,uhatstar,...
                vhstar,phstar,vhatstar,...
                uh,qh,uhat,...
                vh,ph,vhat,...
                GQ1DRef_pts,GQ1DRef_wts,k,k_star,tau)

    Nu = (k+1)*(k+2)/2;
    Nuhat  = k+1;
    dir_vec = GetDirVec(Nuhat); % correct the uhat oritation
    Nustar = (k_star+1)*(k_star+2)/2;
   
    NGQ = length(GQ1DRef_pts);
    num_elements = mymesh.num_elements;
    
    if strcmp(func_type,'1')
         
        % Get Gauss Quadpoints on the square
        [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts);
        % Map Gauss Quadpoints to the reference triangle
        [r_list,s_list] = ABtoRS(a_list,b_list);
        
        if strcmp(pde_ype,'1') % Poission
            
            % Vandermonde matrix for uh, qh, and grad(uh)
            V2D = Vandermonde2D(k,a_list,b_list);
            [V2Dr,V2Ds] = GradVandermonde2D(k,a_list,b_list);
            
            V1D = Vandermonde1D(k,GQ1DRef_pts);% legendre polynomail at GQ points
            
            % Vanermonde matrix for uhstar and qhstar
            VD_star = Vandermonde2D(k_star,a_list,b_list);
            
            
            % compute the estimate terms over all elements
            
            % est1 = (qh* - qh, ph* - ph)
            Est1_elewise_list = zeros(num_elements,1,numeric_t);
            
            % est2 = (q*-qh,ph+grad(vh))
            Est2_elewise_list = zeros(num_elements,1,numeric_t);
            
            % est3 = (qh+grad(uh),p*-ph)
            Est3_elewise_list = zeros(num_elements,1,numeric_t);
            
            % est4 = <qhat-qstar, vh - vhat>
            Est4_elewise_list = zeros(num_elements,1,numeric_t);
            
            %est5 = <uh - uhat, phat-qhstar>
            Est5_elewise_list = zeros(num_elements,1,numeric_t);
            
            extra_term = zeros(num_elements,1,numeric_t);
             
            
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
                
                % uhstar, vhstar at those Gauss points
                uhstar_coeff = uhstar(:,element_idx);
                vhstar_coeff = vhstar(:,element_idx);
                
                uhstar_pts = VD_star*uhstar_coeff;uhstar_pts = reshape(uhstar_pts,[],NGQ);
                
                vhstar_pts = VD_star*vhstar_coeff;vhstar_pts = reshape(vhstar_pts,[],NGQ);
                
                % qhstar, phstar at those Gauss points
                qhstar_coeff_1 = qhstar(1:Nustar,element_idx);
                qhstar_coeff_2 = qhstar(Nustar+1:end,element_idx);

                phstar_coeff_1 = phstar(1:Nustar,element_idx);
                phstar_coeff_2 = phstar(Nustar+1:end,element_idx);

                q1_VD = VD_star*qhstar_coeff_1;q1_VD = reshape(q1_VD,[],NGQ);
                
                q2_VD = VD_star*qhstar_coeff_2;q2_VD = reshape(q2_VD,[],NGQ);
                
                p1_VD = VD_star*phstar_coeff_1;p1_VD = reshape(p1_VD,[],NGQ);
                
                p2_VD = VD_star*phstar_coeff_2;p2_VD = reshape(p2_VD,[],NGQ);
                
                % uh, vh at Gauss points
                
                uh_pts = V2D*uh_coeff ; uh_pts = reshape(uh_pts,[],NGQ);
                vh_pts = V2D*vh_coeff ; vh_pts = reshape(vh_pts,[],NGQ);
                
                
                % qh,ph at Gauss points
                
                qh_pts1 = V2D * (qh_coeff_1);qh_pts1 = reshape(qh_pts1,[],NGQ);
                
                qh_pts2 = V2D * (qh_coeff_2);qh_pts2 = reshape(qh_pts2,[],NGQ);
                
                ph_pts1 = V2D * (ph_coeff_1);ph_pts1 = reshape(ph_pts1,[],NGQ);
                
                ph_pts2 = V2D * (ph_coeff_2);ph_pts2 = reshape(ph_pts2,[],NGQ);
         
            
                % err_1 (q*-qh,p*-ph)
                temp_formula_1 = ((q1_VD-qh_pts1).*(p1_VD-ph_pts1) ...
                                +  (q2_VD-qh_pts2).*(p2_VD-ph_pts2) );
                            
                Est1_elewise_list(element_idx,1) =...
                    Jk*GQ1DRef_wts'*(temp_formula_1.*Jacobian_rs_to_ab )*GQ1DRef_wts;
            
            
                %d/dx = d/dr * dr/dx + d/ds * ds/dx
                grad_vh_1 = V2Dr * vh_coeff * Inv_AffineMap(1,1)...
                            + V2Ds * vh_coeff * Inv_AffineMap(2,1); grad_vh_1 = reshape(grad_vh_1,[],NGQ);
                %d/dy = d/dr * dr/dy + d/ds * ds/dy
                grad_vh_2 = V2Dr * vh_coeff * Inv_AffineMap(1,2)...
                            + V2Ds * vh_coeff * Inv_AffineMap(2,2); grad_vh_2 = reshape(grad_vh_2,[],NGQ);

                grad_uh_1 = V2Dr * uh_coeff * Inv_AffineMap(1,1)...
                            + V2Ds * uh_coeff * Inv_AffineMap(2,1); grad_uh_1 = reshape(grad_uh_1,[],NGQ);

                grad_uh_2 = V2Dr * uh_coeff * Inv_AffineMap(1,2)...
                            + V2Ds * uh_coeff * Inv_AffineMap(2,2); grad_uh_2 = reshape(grad_uh_2,[],NGQ);

                % err_2 = (q*-qh,ph+grad(vh))  

                temp_formula_2 = ...
                     (q1_VD - qh_pts1).*(ph_pts1 + grad_vh_1)...
                    + (q2_VD - qh_pts2).*(ph_pts2 + grad_vh_2);
                
                % err_3 = (qh+grad(uh),p*-ph)    
                temp_formula_3=...
                     (qh_pts1 + grad_uh_1).*(p1_VD-ph_pts1)...
                    +(qh_pts2 + grad_uh_2).*(p2_VD-ph_pts2);

                Est2_elewise_list(element_idx,1) = ...
                    Jk*GQ1DRef_wts'*(temp_formula_2.*Jacobian_rs_to_ab )*GQ1DRef_wts;
                
                Est3_elewise_list(element_idx,1) = ...
                    Jk*GQ1DRef_wts'*(temp_formula_3.*Jacobian_rs_to_ab )*GQ1DRef_wts;
                
                % err_4 <qhat*n - q*n, vh-vhat> 
                % err_5 < uh- uhat, phat*n-p*n>
                ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);
                uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
                [e_list,n1,n2,n3] = GetTriFaceInfo(vertice_list);
                normal_vector = [n1,n2,n3];
                
                % go through all faces
                
                for ii = 1:length(ele_face_idx_list)
                    face_id = ele_face_idx_list(ii);
                    start_id=(face_id-1)*Nuhat+1;
                    end_id = face_id*Nuhat;

                    uhat_coeff = uhat(start_id:end_id,1);
                    vhat_coeff = vhat(start_id:end_id,1);

                    
                    if uhat_dir_list(1,ii) == 0
                        uhat_coeff = uhat_coeff.*dir_vec;
                        vhat_coeff = vhat_coeff.*dir_vec;
                  
                    end
                    uhat_face_pts = V1D * uhat_coeff;
                    vhat_face_pts = V1D * vhat_coeff;

                   
                    [face_r_list,face_s_list] = GetRefFaceQuadPts(ii,GQ1DRef_pts);
                    [face_a_list,face_b_list] = RStoAB(face_r_list,face_s_list);
                    
                    % ustar qhstar 
                    V2Dstar_face = Vandermonde2D(k_star,face_a_list,face_b_list);
                    
                    q_face_pts1 = V2Dstar_face*qhstar_coeff_1;
                    q_face_pts2 = V2Dstar_face*qhstar_coeff_2;
                    p_face_pts1 = V2Dstar_face*phstar_coeff_1;
                    p_face_pts2 = V2Dstar_face*phstar_coeff_2;
                    
                    %%% temp code
                    u_face_pts = V2Dstar_face*uhstar_coeff;
                    
                    % uh, qh
                    V2D_face = Vandermonde2D(k,face_a_list,face_b_list);

                    qh_face_pts1 = V2D_face * qh_coeff_1;
                    qh_face_pts2 = V2D_face * qh_coeff_2; 
                    uh_face_pts  = V2D_face * uh_coeff  ;

                    ph_face_pts1 = V2D_face * ph_coeff_1;
                    ph_face_pts2 = V2D_face * ph_coeff_2;
                    vh_face_pts  = V2D_face * vh_coeff  ;
                    
                    % (qh - qstar)*n
                    qh_q_n = (qh_face_pts1 - q_face_pts1)*normal_vector(1,ii)...
                            +(qh_face_pts2 - q_face_pts2)*normal_vector(2,ii);
                    
                    % uh - uhat    
                    uh_uhat = uh_face_pts - uhat_face_pts;
                    
                    % (ph - pstar)*n
                    ph_p_n = (ph_face_pts1 - p_face_pts1)*normal_vector(1,ii)...
                            +(ph_face_pts2 - p_face_pts2)*normal_vector(2,ii);
                    % vh - vhat 
                    vh_vhat = vh_face_pts - vhat_face_pts;
                    
                    % <qh*n+tau(uh-uhat) - qhstar*n, vh-vhat>
                    temp_formula_4 = ...
                         (qh_q_n + tau*(uh_uhat)).*(vh_vhat);
                         
                    
                    temp_formula_5 = (ph_p_n + tau*(vh_vhat)).*(uh_uhat);

                    Est4_elewise_list(element_idx,1) =  Est4_elewise_list(element_idx,1)+...
                        0.5*e_list(ii)*GQ1DRef_wts'*(temp_formula_4);
                    
                    Est5_elewise_list(element_idx,1) =  Est5_elewise_list(element_idx,1)+...
                        0.5*e_list(ii)*GQ1DRef_wts'*(temp_formula_5);
                    
                    %%%% temp code
                    
                    temp_f_4 = (u_face_pts - uhat_face_pts).*...
                        (p_face_pts1*normal_vector(1,ii)+p_face_pts2*normal_vector(2,ii));
                    
                    extra_term(element_idx,1) =extra_term(element_idx,1)+ 0.5*e_list(ii)*GQ1DRef_wts'*(temp_f_4);
                    

                end
                

                %{
%                 % est4 = (f- Proj_k f,vh* - vh)
%                 L2proj_f_coeff = 1/Jk * Project_F_to_Wh(Jk,vertice_list,...
%                     source_f,k,GQ1DRef_pts,GQ1DRef_wts);
%                 
%                 L2proj_f_pts = V2D*L2proj_f_coeff;  L2proj_f_pts = reshape(L2proj_f_pts,[],NGQ);
%                 
%                 source_f_pts = source_f([x_list,y_list]); source_f_pts = reshape(source_f_pts,[],NGQ);
%                 
%                 temp_formula_4 = (source_f_pts - L2proj_f_pts).*(vhstar_pts - vh_pts );
%                 
%                 Est4_elewise_list(element_idx,1) = ...
%                     Jk*GQ1DRef_wts'*(temp_formula_4.*Jacobian_rs_to_ab )*GQ1DRef_wts;
%                     
%                 % est5 = (uh*-uh, g - Proj_k g)
%                 L2proj_g_coeff = 1/Jk * Project_F_to_Wh(Jk,vertice_list,...
%                     source_g,k,GQ1DRef_pts,GQ1DRef_wts);
%                 
%                 L2proj_g_pts = V2D*L2proj_g_coeff;  L2proj_g_pts = reshape(L2proj_g_pts,[],NGQ);
%                 source_g_pts = source_g([x_list,y_list]); source_g_pts = reshape(source_g_pts,[],NGQ);
%                 
%                 temp_formula_5 = (source_g_pts - L2proj_g_pts).*(uhstar_pts - uh_pts);
%                 
%                 Est5_elewise_list(element_idx,1) = ...
%                     Jk*GQ1DRef_wts'*(temp_formula_5.*Jacobian_rs_to_ab )*GQ1DRef_wts;
%                     
                %}
            end
            
            Est_elewise_list = Est1_elewise_list...
                + Est2_elewise_list +Est3_elewise_list...
                + Est4_elewise_list+Est5_elewise_list+ extra_term;
            
        end
    elseif strcmp(func_type,'2')
        
    else

        error ('This type of error terms calculation has not been implemented yet.')
    end
    
        
end