function [lamh2,lamh_AC,ACh,ACh_Neig_elewise_list]=...
                    EigenvalueFunctional_Maxwell(pde_ype,mymesh,Neig,...
                    wh_Neig,uh_Neig,hat_var_Neig,...
                    GQ1DRef_pts,GQ1DRef_wts,k,tau_t,tau_n,mu,epsilon,omg)
                
       % For Poission Eigen Problem We need the following two parts         
       % part 1: compute lamh = ((wh,curl(uh)) + <what,uhxn> )/(uh,uh)
       % .       namely A =  (wh,curl(uh)) + <what,uhxn>, B = (uh,uh) 
       % part 2: compute adjoint-correction terms
       %         term_1 = 
       %         term_2 = 
       %         term_3 = - < uhxn-uhat_t, what - what >
       % .              = - < uhxn-uhat_t, tau_t*(uhxn-uhat_t) >
                
               
        Nw = (k+1)*(k+2)/2;
        Nu = 2*Nw;
        
        Nuhat = (k+1);
        N_global = 2*Nuhat;
       
        num_elements = mymesh.num_elements;   
        
        Id_mtrix = eye(Nuhat,Nuhat, numeric_t);
        
        
        
        if strcmp(pde_ype,'5') % Maxwell
            B_elewise_list = zeros(num_elements,Neig,numeric_t);
            A_elewise_list = zeros(num_elements,Neig,numeric_t);
            ACh1_elewise_list = zeros(num_elements,Neig,numeric_t);
            ACh2_elewise_list = zeros(num_elements,Neig,numeric_t);
            ACh3_elewise_list = zeros(num_elements,Neig,numeric_t);
            [Aww,Awwr,Awws,Aww3,Bwuhat3] = HDG_LocalMatrix(k,GQ1DRef_pts,GQ1DRef_wts);
            
            for element_idx = 1: num_elements
                %   preparation ----------------------------------------------
                temp_element = mymesh.element_list(element_idx,:);
                vertice_list = mymesh.vertices_list(temp_element(:),:);
                element_faces_list  = mymesh.element_faces_list(element_idx,:);
                Jk = mymesh.Jacobian_list(element_idx);
                uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
                [edge_len_list,~] = GetTriFaceInfo(vertice_list);
                face_type = mymesh.f_type(mymesh.element_faces_list(element_idx,:));

                % local solver matrices
                [Ns,M_Loc] = HDG_MaxwellLocalEquations...
                (Jk,vertice_list,tau_t,tau_n,mu,epsilon,omg,...
                Aww,Awwr,Awws,Aww3,Bwuhat3,uhat_dir_list,face_type);

                
                
                Muw = M_Loc(Nw+1:Nw+Nu,1:Nw);
                Muu = M_Loc(Nw+1:Nw+Nu,Nw+1:Nw+Nu);
                
                
                
                for jj = 1: Neig
                    
                    uh_coeff = uh_Neig(:,element_idx,jj);
                    wh_coeff = wh_Neig(:,element_idx,jj);
                    
                    AA = 0.0;
                    % (div.qh, uh)+<tau*(uh),uh> = -(qh,grad(vh)) + <qh*n+tau(uh),uh>
                    AA = AA + (uh_coeff)'* Muw * wh_coeff; 
                    AA = AA + (uh_coeff)' * Muu * uh_coeff;
                    
                    % <tau_t*(uhxn),uhxn>
                    temp_uh_vh = (uh_coeff)' * Muu * uh_coeff;
                    
                    
                    % <tau*(uh_hat),uhxn>
                    temp_uhat_vh = 0.0 ;
                    % <tau*(uh_hat),uh_hat>
                    temp_uhat_uhat = 0.0;
                    for ii = 1:length(element_faces_list)
                        face_idx = element_faces_list(ii);
                        
                        temp_uhat = hat_var_Neig((face_idx-1)*N_global+1:(face_idx-1)*N_global+Nuhat,jj);

                        Bwuhat = Ns(Nw+1:Nw+Nu,1:Nuhat,ii)  ;
                        
                        temp_1= (uh_coeff)' * Bwuhat * temp_uhat;
                        
                        AA = AA - temp_1;

                        temp_uhat_vh  = temp_uhat_vh + temp_1;
                        
                        

                        Bmu_uhat = tau_t * Id_mtrix * edge_len_list(ii)*0.5;
                        
                        temp_uhat_uhat = temp_uhat_uhat+...
                            temp_uhat' * Bmu_uhat * temp_uhat; % <tau*(uh_hat), uh_hat>

                    end
                    
                    B_elewise_list(element_idx,jj) = Jk*(uh_coeff'*uh_coeff);
                    
                    A_elewise_list(element_idx,jj) = AA;

                    % - <tau*(uh-uh_hat),(uh-uh_hat)>
                    ACh3 = -( (temp_uh_vh  - temp_uhat_vh) - (temp_uhat_vh-temp_uhat_uhat) );
                    
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