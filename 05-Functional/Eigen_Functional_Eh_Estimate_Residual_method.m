function [est_sum,est1,est2,est3,est4]...
        =Eigen_Functional_Eh_Estimate_Residual_method(pde_ype,mymesh,...
                uhstar,qhstar,lamh_tag,...
                uh,qh,uhat,...
                GQ1DRef_pts,GQ1DRef_wts,k,tau)
            
    Nu = (k+1)*(k+2)/2;
    Nuhat  = k+1;
    dir_vec = GetDirVec(Nuhat); % correct the uhat oritation
    
    NGQ = length(GQ1DRef_pts);
    num_elements = mymesh.num_elements;
    
    % Get Gauss Quadpoints on the square
    [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts);
    

    if strcmp(pde_ype,'1') % Poission
        % Vandermonde matrix for uh, qh, and grad(uh)
        V2D = Vandermonde2D(k,a_list,b_list);
        [V2Dr,V2Ds] = GradVandermonde2D(k,a_list,b_list);

        V1D = Vandermonde1D(k,GQ1DRef_pts);% legendre polynomail at GQ points

        % Vanermonde matrix for uhstar and qhstar
        VD_uhstar = Vandermonde2D(k+1,a_list,b_list);
        [V2Dr_star,V2Ds_star] = GradVandermonde2D(k+1,a_list,b_list);
        
        %[VD_qhstar_1,VD_qhstar_2] = RTVandermonde2D(k,a_list,b_list);

        % compute the estimate terms over all elements

        % est1 = (lamhuh* - div.qh, uh* - uh)
        est1 = zeros(num_elements,1,numeric_t);

        % est2 = <(qhat-qh)*n, uh-uh*>
        est2 = zeros(num_elements,1,numeric_t);

        % est3 = (graduh* + qh, graduh* + qh)  --> (q-qh,q-qh)
        est3 = zeros(num_elements,1,numeric_t);
        
        % est6 = lamh(uh* - uh, uh* -uh)
        est4 = zeros(num_elements,1,numeric_t);
        
        uh_L2 = 0.0;

        for element_idx = 1: num_elements
            temp_element = mymesh.element_list(element_idx,:);
            vertice_list = mymesh.vertices_list(temp_element(:),:);
            Jk = mymesh.Jacobian_list(element_idx);
            [~,Inv_AffineMap] = Ref_Tri_Map(Jk,vertice_list);
          

            %% %%%%%%%% compute vhstar, phstar %%%%%%%%%%%%%%%%%%%%%%%

            

            uhstar_coeff = uhstar(:,element_idx);
            uhstar_pts = VD_uhstar*uhstar_coeff;uhstar_pts = reshape(uhstar_pts,[],NGQ);

           
            % graduh*, gradvh*

            grad_uhstar_1 = V2Dr_star * uhstar_coeff * Inv_AffineMap(1,1)...
                        + V2Ds_star * uhstar_coeff * Inv_AffineMap(2,1); grad_uhstar_1 = reshape(grad_uhstar_1,[],NGQ);

            grad_uhstar_2 = V2Dr_star * uhstar_coeff * Inv_AffineMap(1,2)...
                        + V2Ds_star * uhstar_coeff * Inv_AffineMap(2,2); grad_uhstar_2 = reshape(grad_uhstar_2,[],NGQ);


            %% %%%%%%%%%% compute uh,vh, qh,ph related points %%%%%%%%%%

            uh_coeff = uh(:,element_idx);
            
            qh_coeff_1 = qh(1:Nu,element_idx);
            qh_coeff_2 = qh(Nu+1:end,element_idx);
            

            % uh at Gauss points
            uh_pts = V2D*uh_coeff ; uh_pts = reshape(uh_pts,[],NGQ);

            % qh at Gauss points
            qh_pts1 = V2D * (qh_coeff_1);qh_pts1 = reshape(qh_pts1,[],NGQ);
            qh_pts2 = V2D * (qh_coeff_2);qh_pts2 = reshape(qh_pts2,[],NGQ);
            

            % d/dx = d/dr * dr/dx + d/ds * ds/dx
            % d/dy = d/dr * dr/dy + d/ds * ds/dy


            div_qh = V2Dr * qh_coeff_1 * Inv_AffineMap(1,1)+V2Ds * qh_coeff_1 * Inv_AffineMap(2,1)...
                + V2Dr * qh_coeff_2 * Inv_AffineMap(1,2)+ V2Ds * qh_coeff_2 * Inv_AffineMap(2,2);
            div_qh = reshape(div_qh,[],NGQ);


            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % est1 = (lamh*uh* - div.qh, vh* - vh)

            formua_1 = (lamh_tag*uhstar_pts - div_qh ).*(uhstar_pts - uh_pts);
            est1(element_idx,1) =...
                Jk*GQ1DRef_wts'*(formua_1.*Jacobian_rs_to_ab )*GQ1DRef_wts;

            
            % est 5 = (graduh* + qh, graduh*+qh)

            formula_5 =  (grad_uhstar_1+qh_pts1).*(grad_uhstar_1+qh_pts1) + (grad_uhstar_2+qh_pts2).*(grad_uhstar_2+qh_pts2) ;
            est3(element_idx,1) =...
                Jk*GQ1DRef_wts'*(formula_5.*Jacobian_rs_to_ab )*GQ1DRef_wts;

            formula_6 = lamh_tag*(uhstar_pts - uh_pts).^2;
            est4(element_idx,1)=...
                Jk*GQ1DRef_wts'*(formula_6.*Jacobian_rs_to_ab )*GQ1DRef_wts;
            
            formula_7 = uh_pts.^2;
            
            uh_L2 = uh_L2+ Jk*GQ1DRef_wts'*(formula_7.*Jacobian_rs_to_ab )*GQ1DRef_wts;

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
                

                if uhat_dir_list(1,ii) == 0
                    uhat_coeff = uhat_coeff.*dir_vec;
                    
                end

                uhat_face_pts = V1D * uhat_coeff;
                
                [face_r_list,face_s_list] = GetRefFaceQuadPts(ii,GQ1DRef_pts);
                [face_a_list,face_b_list] = RStoAB(face_r_list,face_s_list);


                V2Dstar_face = Vandermonde2D(k+1,face_a_list,face_b_list);

                
                uhstar_face_pts = V2Dstar_face*uhstar_coeff;

                V2D_face = Vandermonde2D(k,face_a_list,face_b_list);

                uh_face_pts  = V2D_face * uh_coeff  ;

                % uh - uhat    
                uh_uhat = uh_face_pts - uhat_face_pts;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % est3 = <(qhat-qh)*n, vh-vh*>

                formula_3 = tau*(uh_uhat).*(uh_face_pts - uhstar_face_pts );

                est2(element_idx,1) =  est2(element_idx,1)+...
                    0.5*e_list(ii)*GQ1DRef_wts'*(formula_3);

            end

        end

        est_sum = (est4 + est3 - 2.0*(est2 + est1))/uh_L2;

    end


    
    
            
            
    
end