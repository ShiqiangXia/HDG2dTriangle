function [lamh2,lamh_AC,ACh,ACh_Neig_elewise_list]=...
                    EigenvalueFunctional_Elliptic(pde_ype,mymesh,Neig,...
                    uh_Neig,qh_Neig,uhat_Neig,...
                    GQ1DRef_pts,GQ1DRef_wts,k,tau)
                
       % For Poission Eigen Problem We need the following two parts         
       % part 1: compute lamh = ( -(qh,grad(uh)) + <uh, qhat*n> )/(uh,uh)
       % .       namely A =  -(qh,grad(uh)) + <uh, qhat*n>, B = (uh,uh) 
       % part 2: compute adjoint-correction terms
       %         term_1 = - ( (qh,qh+grad(uh)) + <uhat -uh,qh*n> )
       %         term_2 = - < uhat, qhat*n >
       %         term_3 = - < uh-uhat, qhat*n - qh*n >
       % .              = - < uh-uhat, tau*(uh-uhat) >
                
                
        Nu = (k+1)*(k+2)/2;
        Nq = 2*Nu;
        Nuhat = k+1;
       
        num_elements = mymesh.num_elements;   
        
        Id_mtrix = eye(Nuhat,Nuhat, numeric_t);
        
        
        
        if strcmp(pde_ype,'1') % Poission
            B_elewise_list = zeros(num_elements,Neig,numeric_t);
            A_elewise_list = zeros(num_elements,Neig,numeric_t);
            ACh1_elewise_list = zeros(num_elements,Neig,numeric_t);
            ACh2_elewise_list = zeros(num_elements,Neig,numeric_t);
            ACh3_elewise_list = zeros(num_elements,Neig,numeric_t);
            [Aqq,Auur,Auus,Auu3,Buuhat3] = HDG_LocalMatrix(k,GQ1DRef_pts,GQ1DRef_wts);
            
            for element_idx = 1: num_elements
                %   preparation ----------------------------------------------
                temp_element = mymesh.element_list(element_idx,:);
                vertice_list = mymesh.vertices_list(temp_element(:),:);
                element_faces_list  = mymesh.element_faces_list(element_idx,:);
                Jk = mymesh.Jacobian_list(element_idx);
                uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
                [edge_len_list,~] = GetTriFaceInfo(vertice_list);

                % local solver matrices
                [Ns,M_Loc] = HDG_PoissionLocalEquations...
                (Jk,vertice_list,tau,Aqq,Auur,Auus,Auu3,Buuhat3,uhat_dir_list);

                %  (div.qh, vh) + <tau*(uh-uh_hat),vh> = (div.qh, uh)+<qhat*n,uh>
                
                Mwq = M_Loc(Nq+1:end,1:Nq);
                Mwu = M_Loc(Nq+1:end,Nq+1:end);
                
                % -(qh,qh) + (uh, div.qh)+<uhat, qh> = -( qh,qh+grad(uh) ) - <uhat -uh,qh*n>
                Mrq = M_Loc(1:Nq,1:Nq);
                Mru = - Mwq';
                
                for jj = 1: Neig
                    
                    uh_coeff = uh_Neig(:,element_idx,jj);
                    qh_coeff = qh_Neig(:,element_idx,jj);
                    
                    AA = 0.0;
                    % (div.qh, uh)+<tau*(uh),uh> = -(qh,grad(vh)) + <qh*n+tau(uh),uh>
                    AA = AA + (uh_coeff)'* Mwq * qh_coeff; 
                    AA = AA + (uh_coeff)' * Mwu * uh_coeff;
                    
                    % <tau*(uh),uh>
                    temp_uh_vh = (uh_coeff)' * Mwu * uh_coeff;
                    
                    ACh1 = 0.0;
                    %  - (qh,qh)
                    ACh1 = ACh1 - qh_coeff' * Mrq * qh_coeff;
                    % - (-(uh, div.qh))
                    ACh1 = ACh1 - qh_coeff' * Mru * uh_coeff;
                    
                    % <tau*(uh_hat),uh>
                    temp_uhat_vh = 0.0 ;
                    
                    % <qh*n + tau*(uh-uh_hat), vh_hat>
                    ACh2 = 0.0;
                    temp_uhuhat_vhat = 0.0;
                    for ii = 1:length(element_faces_list)
                        face_idx = element_faces_list(ii);
                        temp_uhat = uhat_Neig((face_idx-1)*Nuhat+1:face_idx*Nuhat,jj);

                        Bwuhat = Ns(Nq+1:end,:,ii)  ;
                        temp_1= (uh_coeff)' * Bwuhat * temp_uhat;
                        
                        AA = AA - temp_1;

                        temp_uhat_vh  = temp_uhat_vh + temp_1;
                        
                        Bruhat = -Ns(1:Nq,:,ii)  ;
                        %  - <uhat, qh>
                        temp_2 = qh_coeff' * Bruhat * temp_uhat;
                        ACh1 = ACh1 - temp_2;
                        
                        % <qh*n + tau*(uh), vh_hat>  - tau<uh_hat,vh_hat>
                        %Bmu_q = (-Ns(1:Nq,:,ii))';
                        %Bmu_u = Ns(Nq+1:end,:,ii)';

                        Bmu_uhat = tau * Id_mtrix * edge_len_list(ii)*0.5;
                        temp_3 = temp_2;
                        %temp_3 = temp_uhat' *  Bmu_q * qh_coeff; % <qh*n, uh_hat>
                        temp_5 = temp_1;
                        %temp_5 = temp_uhat' * Bmu_u * uh_coeff; % <tau*(uh), uh_hat>
                        temp_4 = temp_uhat' * Bmu_uhat * temp_uhat; % <tau*(uh_hat), uh_hat>
                        
                        %<qh*n + tau*(uh-uhat), uh_hat>
                        ACh2 = ACh2 + temp_3 + temp_5 - temp_4;
                        
                        temp_uhuhat_vhat = temp_uhuhat_vhat +temp_5 - temp_4;

                        
                    end
                    
                    B_elewise_list(element_idx,jj) = Jk*(uh_coeff'*uh_coeff);
                    
                    A_elewise_list(element_idx,jj) = AA;
                    
                    ACh1_elewise_list(element_idx,jj) =  ACh1;
                    
                    ACh2_elewise_list(element_idx,jj) = -ACh2;
                    
                    % - <tau*(uh-uh_hat),(uh-uh_hat)>
                    ACh3 = +( (temp_uh_vh  - temp_uhat_vh) - temp_uhuhat_vhat );
                    
                    ACh3_elewise_list(element_idx,jj) = ACh3;
                    
                end
                %----------------------------------------------------------

                
            end
            
            
            
            B_Neig = sum(B_elewise_list);
            A_Neig = sum(A_elewise_list);
            
            ACh_up    = sum(ACh1_elewise_list+ACh2_elewise_list+ACh3_elewise_list);
            
            ACh_Neig_elewise_list = ACh3_elewise_list;
            
            
        end
        
        
        lamh2 = A_Neig./B_Neig;
        ACh = ACh_up./B_Neig;
        lamh_AC = lamh2  + ACh;
        
        

                
               
end