function [est_sum,est1,est2,est3,est4]=Functional_Eh_Estimate_Residual_method(func_type,pde_ype,mymesh,...
                source_f,...
                vhstar,phstar,...
                uh,qh,uhat,...
                vh,ph,vhat,...
                GQ1DRef_pts,GQ1DRef_wts,k,tau)
            
    Nu = (k+1)*(k+2)/2;
    Nuhat  = k+1;
    dir_vec = GetDirVec(Nuhat); % correct the uhat oritation
    
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
            VD_uhstar = Vandermonde2D(k+1,a_list,b_list);
            [V2Dr_star,V2Ds_star] = GradVandermonde2D(k+1,a_list,b_list);
            [VD_qhstar_1,VD_qhstar_2] = RTVandermonde2D(k,a_list,b_list);
            
            % compute the estimate terms over all elements
            
            % est1 = (f- div.qh, vh* - vh)
            est1 = zeros(num_elements,1,numeric_t);
            
            % est2 = (qh+grad uh, ph* - ph)
            est2 = zeros(num_elements,1,numeric_t);
            
            % est3 = <(qhat-qh)*n, vh-vh*>
            est3 = zeros(num_elements,1,numeric_t);
            
            % est4 = <uh - uhat, (phat + grad uh*)*n>
            est4 = zeros(num_elements,1,numeric_t);
            
            for element_idx = 1: num_elements
                temp_element = mymesh.element_list(element_idx,:);
                vertice_list = mymesh.vertices_list(temp_element(:),:);
                Jk = mymesh.Jacobian_list(element_idx);
                [~,Inv_AffineMap] = Ref_Tri_Map(Jk,vertice_list);
                % Gauss points on any element
                [x_list,y_list] = RStoXY(r_list,s_list,Jk,vertice_list);
                
                %% %%%%%%%% compute vhstar, phstar %%%%%%%%%%%%%%%%%%%%%%%
                
                % vhstar at those Gauss points
                vhstar_coeff = vhstar(:,element_idx);
                vhstar_pts = VD_uhstar*vhstar_coeff;vhstar_pts = reshape(vhstar_pts,[],NGQ);
                
                % phstar at those Gauss points
                phstar_coeff = phstar(:,element_idx);
                p1_VD = VD_qhstar_1*phstar_coeff;p1_VD = reshape(p1_VD,[],NGQ);
                p2_VD = VD_qhstar_2*phstar_coeff;p2_VD = reshape(p2_VD,[],NGQ);
                
                %% %%%%%%%%%% compute uh,vh, qh,ph related points %%%%%%%%%%
                
                uh_coeff = uh(:,element_idx);
                vh_coeff = vh(:,element_idx);
                qh_coeff_1 = qh(1:Nu,element_idx);
                qh_coeff_2 = qh(Nu+1:end,element_idx);
                ph_coeff_1 = ph(1:Nu,element_idx);
                ph_coeff_2 = ph(Nu+1:end,element_idx);
                
                % vh at Gauss points
                vh_pts = V2D*vh_coeff ; vh_pts = reshape(vh_pts,[],NGQ);
                
                % qh,ph at Gauss points
                qh_pts1 = V2D * (qh_coeff_1);qh_pts1 = reshape(qh_pts1,[],NGQ);
                qh_pts2 = V2D * (qh_coeff_2);qh_pts2 = reshape(qh_pts2,[],NGQ);
                ph_pts1 = V2D * (ph_coeff_1);ph_pts1 = reshape(ph_pts1,[],NGQ);
                ph_pts2 = V2D * (ph_coeff_2);ph_pts2 = reshape(ph_pts2,[],NGQ);
                
                % d/dx = d/dr * dr/dx + d/ds * ds/dx
                % d/dy = d/dr * dr/dy + d/ds * ds/dy
                
                grad_uh_1 = V2Dr * uh_coeff * Inv_AffineMap(1,1)...
                            + V2Ds * uh_coeff * Inv_AffineMap(2,1); grad_uh_1 = reshape(grad_uh_1,[],NGQ);

                grad_uh_2 = V2Dr * uh_coeff * Inv_AffineMap(1,2)...
                            + V2Ds * uh_coeff * Inv_AffineMap(2,2); grad_uh_2 = reshape(grad_uh_2,[],NGQ);
                        
                div_qh = V2Dr * qh_coeff_1 * Inv_AffineMap(1,1)+V2Ds * qh_coeff_1 * Inv_AffineMap(2,1)...
                    + V2Dr * qh_coeff_2 * Inv_AffineMap(1,2)+ V2Ds * qh_coeff_2 * Inv_AffineMap(2,2);
                div_qh = reshape(div_qh,[],NGQ);
                
                source_f_pts = source_f([x_list,y_list]); source_f_pts = reshape(source_f_pts,[],NGQ);
                
                
                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                formua_1 = (source_f_pts - div_qh ).*(vhstar_pts - vh_pts);
                est1(element_idx,1) =...
                    Jk*GQ1DRef_wts'*(formua_1.*Jacobian_rs_to_ab )*GQ1DRef_wts;
                
                formua_2 = (qh_pts1 + grad_uh_1).*(p1_VD - ph_pts1) ...
                    + (qh_pts2 + grad_uh_2).*(p2_VD -ph_pts2 );
                
                est2(element_idx,1) =...
                    Jk*GQ1DRef_wts'*(formua_2.*Jacobian_rs_to_ab )*GQ1DRef_wts;
                
                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                    
                    
                    V2Dstar_face = Vandermonde2D(k+1,face_a_list,face_b_list);
                    [V2Dr_star_face,V2Ds_star_face] = GradVandermonde2D(k+1,face_a_list,face_b_list);
                    
                    % phstar = - grad vh*
                    
                    phstar_face_1 =  - (V2Dr_star_face * vhstar_coeff * Inv_AffineMap(1,1)...
                        + V2Ds_star_face * vhstar_coeff * Inv_AffineMap(2,1));phstar_face_1 = reshape(phstar_face_1,[],NGQ);
                    
                    phstar_face_2 = - (V2Dr_star_face * vhstar_coeff * Inv_AffineMap(1,2)...
                        +V2Ds_star_face * vhstar_coeff * Inv_AffineMap(2,2)); phstar_face_2 = reshape(phstar_face_2,[],NGQ);
                    
                    
                    vhstar_face_pts = V2Dstar_face*vhstar_coeff;
                    
                    V2D_face = Vandermonde2D(k,face_a_list,face_b_list);
                    qh_face_pts1 = V2D_face * qh_coeff_1;
                    qh_face_pts2 = V2D_face * qh_coeff_2; 
                    uh_face_pts  = V2D_face * uh_coeff  ;

                    ph_face_pts1 = V2D_face * ph_coeff_1;
                    ph_face_pts2 = V2D_face * ph_coeff_2;
                    vh_face_pts  = V2D_face * vh_coeff  ;
                    
                    % uh - uhat    
                    uh_uhat = uh_face_pts - uhat_face_pts;
                    % vh - vhat 
                    vh_vhat = vh_face_pts - vhat_face_pts;
                    
                    ph_pstar_pts = (ph_face_pts1- phstar_face_1)*normal_vector(1,ii)...
                        +(ph_face_pts2-phstar_face_2)*normal_vector(2,ii);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % est3 = <(qhat-qh)*n, vh-vh*>
                     
                    formula_3 = tau*(uh_uhat).*(vh_face_pts - vhstar_face_pts );
                    
                    est3(element_idx,1) =  est3(element_idx,1)+...
                        0.5*e_list(ii)*GQ1DRef_wts'*(formula_3);
                    
                    % est4 = <uh - uhat, (phat + grad uh*)*n>
                    formula_4 = uh_uhat.*(ph_pstar_pts+tau*vh_vhat);
                    
                    est4(element_idx,1) =  est4(element_idx,1)+...
                        0.5*e_list(ii)*GQ1DRef_wts'*(formula_4);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                end

            end
            
            est_sum = est1 +est2 +est3 +est4;
            
        end
        
        
    else
    end
    
    
            
            
    
end