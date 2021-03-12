function [Jh,Jh_AC,ACh,ACh_elewise_list] = LinearFunctional(func_type,pde_ype,mymesh,...
        uh,qh,uhat,source_f,...
        vh,ph,vhat,source_g,...
        GQ1DRef_pts,GQ1DRef_wts,k,tau,post_flag)
    % compute the linear functiona and adjoint correction terms
    
    % notice that many integral terms have the form 
    % (a_h, b_h) = coefficient_a^T  * Mat_{a,b} * coefficient_b
    % coefficient_a, coefficient_b are the DG solutions (puts) 
    %  Mat_{a,b} are matrices we have already computed when solving the PDE.
    % Therefore, to compute these terms, we just need to call the local solver
    % equations to obtain these matrices and do matrix multiplication. 
    
    Nu = (k+1)*(k+2)/2;
    Nq = 2*Nu;
    Nuhat = k+1;
    NGQ = length(GQ1DRef_pts);
    
    num_elements = mymesh.num_elements; 
    
    Id_mtrix = eye(Nuhat,Nuhat, numeric_t);
    
    Jh = 0.0;
    
    % Get Gauss Quadpoints on the square
    [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts);
    % Map Gauss Quadpoints to the reference triangle
    [r_list,s_list] = ABtoRS(a_list,b_list);
    
    if post_flag == 0
        V2D = Vandermonde2D(k,a_list,b_list);
    elseif post_flag == 1
        error('Post-processed has not implemented yet');
    end
    
    
    if strcmp(func_type,'1')
     % functional J(u) = (u,g)
     % Jh(uh) = (uh,g);

     if strcmp(pde_ype,'1') % Poission
         
         ACh1_elewise_list = zeros(num_elements,1,numeric_t);
         ACh2_elewise_list = zeros(num_elements,1,numeric_t);
         ACh3_elewise_list = zeros(num_elements,1,numeric_t);
         ACh4_elewise_list = zeros(num_elements,1,numeric_t);
         [Aqq,Auur,Auus,Auu3,Buuhat3] = HDG_LocalMatrix(k,GQ1DRef_pts,GQ1DRef_wts);
         
         for element_idx = 1: num_elements
             
             %   preparation ----------------------------------------------
             temp_element = mymesh.element_list(element_idx,:);
             vertice_list = mymesh.vertices_list(temp_element(:),:);
             ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);
             Jk = mymesh.Jacobian_list(element_idx);
             uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
             [edge_len_list,~] = GetTriFaceInfo(vertice_list);
             
             uh_coeff = uh(:,element_idx);
             vh_coeff = vh(:,element_idx);
             qh_coeff = qh(:,element_idx);
             ph_coeff = ph(:,element_idx);
             
             % local solver matrices
             [Ns,M_Loc] = HDG_PoissionLocalEquations...
                 (Jk,vertice_list,tau,Aqq,Auur,Auus,Auu3,Buuhat3,uhat_dir_list);
 
             
             % Jh(uh) = (g,uh) ----------------------------------------------
             % Gauss points on any element
             [x_list,y_list] = RStoXY(r_list,s_list,Jk,vertice_list);
             
             g_VD = source_g([x_list,y_list]);
             
             uh_pts = V2D * (uh_coeff); % uh on Gauss points
             uh_pts = reshape(uh_pts,[],NGQ);
             
             Jh = Jh + Jk*GQ1DRef_wts'*(g_VD.*uh_pts.*Jacobian_rs_to_ab)*GQ1DRef_wts;
             
             % adjoint correction terms ----------------------------------------------
            
             % ACh1 = (f,vh) -(div.qh, vh) - <tau*(uh-uh_hat),vh>
             ACh1 = 0.0;
             
             f_VD = source_f([x_list,y_list]);
             vh_pts = V2D * (vh_coeff); % uh on Gauss points
             vh_pts = reshape(vh_pts,[],NGQ);
             
             
             
                %(f,vh)
            
             ACh1 = ACh1 + Jk*GQ1DRef_wts'*(f_VD.*vh_pts.*Jacobian_rs_to_ab)*GQ1DRef_wts;
             
                % - (div.qh, vh) - <tau*(uh),vh>
             Mwq = M_Loc(Nq+1:end,1:Nq);
             Mwu = M_Loc(Nq+1:end,Nq+1:end);
             ACh1 = ACh1 - (vh_coeff)'* Mwq * qh_coeff;
             ACh1 = ACh1 - (vh_coeff)' * Mwu * uh_coeff; 
             
             temp_uh_vh = (vh_coeff)' * Mwu * uh_coeff;
             
                % + <tau*(uh_hat),vh>
                
             temp_uhat_vh = 0.0 ;
             for ii = 1:length(ele_face_idx_list)
                 face_idx = element_faces_list(ii);
                 temp_uhat = uhat((face_idx-1)*Nuhat+1:face_idx*Nuhat);
                 
                 Bwuhat = Ns(Nq+1:end,:,ii)  ;
                 temp_1= (vh_coeff)' * Bwuhat * temp_uhat;
                 ACh1 = ACh1 + temp;
                 
                 temp_uhat_vh  = temp_uhat_vh + temp_1;
             end
             
             ACh1_elewise_list(element_idx,1) = ACh1;
             
             % ACh2 = (qh,ph) - (uh, div.ph)+<uhat, ph>
             
             ACh2 = 0.0;
             
             Mrq = M_Loc(1:Nq,1:Nq);
             Mru = - Mwq';
             ACh2 = ACh2 + ph_coeff' * Mrq * qh_coeff;
             ACh2 = ACh2 + ph_coeff' * Mru * uh_coeff;
             
             for ii = 1:length(ele_face_idx_list)
                 face_idx = element_faces_list(ii);
                 temp_uhat = uhat((face_idx-1)*Nuhat+1:face_idx*Nuhat);
                 Bruhat = -Ns(1:Nq,:,ii)  ;
                 ACh2 = ACh2 + ph_coeff' * Bruhat * temp_uhat;
             end
             
             ACh2_elewise_list(element_idx,1) = ACh2;
             
             
             % ACh3 = <qh*n + tau*(uh-uh_hat), vh_hat>
             
             ACh3 = 0.0;
             temp_uhuhat_vhat = 0.0;
             for ii = 1:length(ele_face_idx_list)
                 face_idx = element_faces_list(ii);
                 
                 temp_vhat = vhat((face_idx-1)*Nuhat+1:face_idx*Nuhat);
                 temp_uhat = uhat((face_idx-1)*Nuhat+1:face_idx*Nuhat);
                 
                 % <qh*n + tau*(uh), vh_hat>  - tau<uh_hat,vh_hat>
                 Bmu_q = (-Ns(1:Nq,:,ii))';
                 Bmu_u = Ns(Nq+1:end,:,ii)';
                 
                 Bmu_uhat = tau * Id_mtrix * edge_len_list(ii)*0.5;
                 
                 temp_3 = temp_vhat' *  Bmu_q * qh_coeff; % <qh*n, vh_hat>
                 temp_5 = temp_vhat' * Bmu_u * uh_coeff; % <tau*(uh), vh_hat>
                 temp_4 = temp_vhat' * Bmu_uhat * temp_uhat; % <tau*(uh_hat), vh_hat>
                 
                 ACh3 = ACh3 + temp_3 + temp_5 - temp_4;
                 
                 temp_uhuhat_vhat = temp_uhuhat_vhat +temp_5 - temp_4;
                 
             end
             
             ACh3_elewise_list(element_idx,1) = ACh3;
             

             
             % ACH4 = - <tau*(uh-uh_hat),(vh-vh_hat)>
             
             ACh4 = -( (temp_uh_vh  - temp_uhat_vh) - temp_uhuhat_vhat );
             
             ACh4_elewise_list(element_idx,1) = ACh4  ;           
             
         end
         
         ACh_elewise_list = ACh1_elewise_list +  ACh2_elewise_list...
                            + ACh3_elewise_list + ACh4_elewise_list;
         
     end

    else
    end
    
    
    ACh = sum(ACh_elewise_list);
    Jh_AC = Jh + ACh;
    
end