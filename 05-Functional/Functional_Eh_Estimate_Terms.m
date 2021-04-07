function [Est_elewise_list,Est1_elewise_list,...
        Est2_elewise_list,Est3_elewise_list,...
        Est4_elewise_list,Est5_elewise_list] ...
                = Functional_Eh_Estimate_Terms(func_type,pde_ype,mymesh,...
                uhstar,qhstar,source_f,...
                vhstar,phstar,source_g,...
                uh,qh,uhat,...
                vh,ph,vhat,...
                GQ1DRef_pts,GQ1DRef_wts,k)

    Nu = (k+1)*(k+2)/2;
   
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
            
            % Vanermonde matrix for uhstar and qhstar
            VD_uhstar = Vandermonde2D(k+1,a_list,b_list);
            [VD_qhstar_1,VD_qhstar_2] = RTVandermonde2D(k,a_list,b_list);
            
            % compute the estimate terms over all elements
            
            % est1 = (qh* - qh, ph* - ph)
            Est1_elewise_list = zeros(num_elements,1,numeric_t);
            
            % est2 = (q*-qh,ph+grad(vh))
            Est2_elewise_list = zeros(num_elements,1,numeric_t);
            
            % est3 = (qh+grad(uh),p*-ph)
            Est3_elewise_list = zeros(num_elements,1,numeric_t);
            
            % est4 = (f- Proj_k f,vh* - vh)
            Est4_elewise_list = zeros(num_elements,1,numeric_t);
            
            %est5 = (uh*-uh, g - Proj_k g)
            Est5_elewise_list = zeros(num_elements,1,numeric_t);
             
            
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
                
                uhstar_pts = VD_uhstar*uhstar_coeff;uhstar_pts = reshape(uhstar_pts,[],NGQ);
                
                vhstar_pts = VD_uhstar*vhstar_coeff;vhstar_pts = reshape(vhstar_pts,[],NGQ);
                
                % qhstar, phstar at those Gauss points
                qhstar_coeff = qhstar(:,element_idx);

                phstar_coeff = phstar(:,element_idx);

                q1_VD = VD_qhstar_1*qhstar_coeff;q1_VD = reshape(q1_VD,[],NGQ);
                
                q2_VD = VD_qhstar_2*qhstar_coeff;q2_VD = reshape(q2_VD,[],NGQ);
                
                p1_VD = VD_qhstar_1*phstar_coeff;p1_VD = reshape(p1_VD,[],NGQ);
                
                p2_VD = VD_qhstar_2*phstar_coeff;p2_VD = reshape(p2_VD,[],NGQ);
                
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

                
                % est4 = (f- Proj_k f,vh* - vh)
                L2proj_f_coeff = 1/Jk * Project_F_to_Wh(Jk,vertice_list,...
                    source_f,k,GQ1DRef_pts,GQ1DRef_wts);
                
                L2proj_f_pts = V2D*L2proj_f_coeff;  L2proj_f_pts = reshape(L2proj_f_pts,[],NGQ);
                
                source_f_pts = source_f([x_list,y_list]); source_f_pts = reshape(source_f_pts,[],NGQ);
                
                temp_formula_4 = (source_f_pts - L2proj_f_pts).*(uhstar_pts - uh_pts );
                
                Est4_elewise_list(element_idx,1) = ...
                    Jk*GQ1DRef_wts'*(temp_formula_4.*Jacobian_rs_to_ab )*GQ1DRef_wts;
                    
                % est5 = (uh*-uh, g - Proj_k g)
                L2proj_g_coeff = 1/Jk * Project_F_to_Wh(Jk,vertice_list,...
                    source_g,k,GQ1DRef_pts,GQ1DRef_wts);
                
                L2proj_g_pts = V2D*L2proj_g_coeff;  L2proj_g_pts = reshape(L2proj_g_pts,[],NGQ);
                source_g_pts = source_g([x_list,y_list]); source_g_pts = reshape(source_g_pts,[],NGQ);
                
                temp_formula_5 = (source_g_pts - L2proj_g_pts).*(vhstar_pts - vh_pts);
                
                Est5_elewise_list(element_idx,1) = ...
                    Jk*GQ1DRef_wts'*(temp_formula_5.*Jacobian_rs_to_ab )*GQ1DRef_wts;
                    
                
            end
            
            Est_elewise_list = Est1_elewise_list + Est2_elewise_list +Est3_elewise_list + Est4_elewise_list+Est5_elewise_list;
            
        end
    elseif strcmp(func_type,'2')
        
    else

        error ('This type of error terms calculation has not been implemented yet.')
    end
    
        
end