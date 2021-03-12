function [uh,qh,uh_hat] = HDG_GlobalSolver(pb, mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,tau,source_f,uD,uN)
    
    % this global solver works for any source problem with dirichlet boundary
    % namely solve:
    % sum <qh*n+tau(uh-uhat),mu>_pK = 0
    %
    % Neuman boundary is not implemented yet. 
    %
    % for different type of PDE,the difference is to call different local solver 
    % matrices in Step 1 and Step2. 
    % the rest assmeble global matrix, solve and recover part is the same!
    
    %% ----- Step 0 Set up parameters ---------------------------------------
    Nu = (k+1)*(k+2)/2;
    Nq = 2*Nu;
    Nuhat = k+1;
    
    num_elements = mymesh.num_elements;
    num_faces = mymesh.num_faces;
    
    %% ----- step 1 get Local Matrices on ref element ----------------------------------
    
    if strcmp(pb,'1')
        % Poission Problem
        [Aqq,Auur,Auus,Auu3,Buuhat3] = HDG_LocalMatrix(k,GQ1DRef_pts,GQ1DRef_wts);
    else
        error('HDG method for problem type %s has not implemented yet ', pb);
    end
    
    %% ----- step 2 Get Local Solver matrices -------------------------
    
    %%% Local equation
    %
    %          qh
    %  M_loc* (  ) = Sum: N_i * uhat_Fi +Proj_f
    %          uh
    % >>
    %  (qh = (Q   * uhat   + (Qw
    %   uh)   U)              Uw)
    %  
    %%%%
    List_LocSol_f = zeros(Nq+Nu,num_elements,numeric_t); % Local solver Qw,Uw  <----f
    List_LocSol = zeros(Nq+Nu,Nuhat,num_elements,3,numeric_t); % Local solver Q, U <----uhat
    List_Ns = zeros(Nq+Nu,Nuhat,num_elements,3,numeric_t); % helper for Global Eqs
    % go through all the elements to build local solver matrices
    
    for element_idx = 1: num_elements
        
        temp_element = mymesh.element_list(element_idx,:);
        
        vertice_list = mymesh.vertices_list(temp_element(:),:);
        Jk = mymesh.Jacobian_list(element_idx);
        
        uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
        
        % Local equation matrices
        if strcmp(pb,'1')
            [Ns,M_Loc] = HDG_PoissionLocalEquations(...
                                    Jk,vertice_list,tau,...
                                       Aqq,Auur,Auus,Auu3,Buuhat3,uhat_dir_list);
        end
        
        % Local solver Q,U
        % M_loc^-1 * N_i
        for jj = 1:3
            List_LocSol(:,:,element_idx,jj) = M_Loc\Ns(:,:,jj);
            List_Ns(:,:,element_idx,jj) = [-Ns(1:Nq,:,jj);Ns(Nq+1:end,:,jj)];
        end
        
%         List_LocSol(:,:,element_idx,2) = M_Loc\N_2;
%         List_Ns(:,:,element_idx,2) = [-N_2(1:Nq,:);N_2(Nq+1:end,:)];
%         
%         List_LocSol(:,:,element_idx,3) = M_Loc\N_3;
%         List_Ns(:,:,element_idx,3) = [-N_3(1:Nq,:);N_3(Nq+1:end,:)];
%         
        % Compute the projection of source_f
        Proj_f = zeros(Nq+Nu,1,numeric_t);
        
        Proj_f(Nq+1:end) = Project_F_to_Wh(Jk,vertice_list,...
            source_f,k,...
            GQ1DRef_pts,GQ1DRef_wts);
        
        % Local solver Qw * f , Uw * f
        % M_loc^-1 * Proj_f
        List_LocSol_f(:,element_idx) = M_Loc\Proj_f;

    end
    
    %% ----- step 3 Assemble Global Matrix and right hand side --------
    
    % Preparation
    Global_M = numeric_t(sparse(Nuhat*num_faces,Nuhat*num_faces));
    Global_b = numeric_t(sparse(Nuhat*num_faces,1));
    Id_mtrix = eye(Nuhat,Nuhat, numeric_t);
    
    for element_idx = 1: num_elements
        
        ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);
        temp_element = mymesh.element_list(element_idx,:);
        vertice_list = mymesh.vertices_list(temp_element(:),:);
        
        uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
        
        [edge_len_list,~] = GetTriFaceInfo(vertice_list);
        for ii = 1:length(ele_face_idx_list)
            face_id = ele_face_idx_list(ii);
            
            start_id=(face_id-1)*Nuhat+1;
            end_id = face_id*Nuhat;
         
            bdry_flag = mymesh.f_type(face_id);
            
            if bdry_flag == 0 % interior face
                
                Bd_Int_mat = List_Ns(:,:,element_idx,ii)';
                
                % put matrix block at the right position
                %<qh*n+tau*uh,mu>
                for jj = 1:length(ele_face_idx_list)
                    temp_id = ele_face_idx_list(jj);
                    temp_start = (temp_id-1)*Nuhat+1;
                    temp_end = temp_id*Nuhat;
                    
                    Global_M(start_id:end_id,temp_start:temp_end) =...
                        Global_M(start_id:end_id,temp_start:temp_end)+...
                       Bd_Int_mat * List_LocSol(:,:,element_idx,jj) ;
                end
                
                Global_b(start_id:end_id,1) =Global_b(start_id:end_id,1)-Bd_Int_mat*List_LocSol_f(:,element_idx);
                
                % -<tau uhat,mu>
                Global_M(start_id:end_id,start_id:end_id) = ...
                    Global_M(start_id:end_id,start_id:end_id) - tau*Id_mtrix*edge_len_list(ii)*0.5;
                     
            elseif bdry_flag == 1 % dirichlet boundary 
                % dirichlet boundary face 
                % just simply enforce boundary conditions
                
                Jk = mymesh.Jacobian_list(element_idx);

                Global_M(start_id:end_id,start_id:end_id)= Id_mtrix*edge_len_list(ii)*0.5;
                
                Global_b(start_id:end_id,1) = Project_F_to_Face(Jk,vertice_list,...
                    ii,uhat_dir_list(1,ii),edge_len_list(ii),...
                    k,uD,GQ1DRef_pts,GQ1DRef_wts);
                
                
            elseif bdry_flag == 2 % neuman boundary
                % Will implement later
                error(' Boundary type not implemented yet.')
            else 
                error('Wrong boundary type.')
            end
            
        end
    end
    
    %% ------ Step 4: Solve the Global Equals ------------------------
    
    uh_hat = Global_M\Global_b;
    
    %% ------ Step 5: Recover Local Solutions --------------------------
    
    %%% Local solvers
    %  ( qh )  =   Q*uhat + Qw*f
    %  ( uh )  =   U*uhat + Uw*f
    %
    % List_LocSol --> (Q,U)
    % List_LocSol_f-->(Qw*f,Uw*f)
    %%%%
    
    Result_matrix = zeros(Nq+Nu,num_elements, numeric_t);
    
    for ele_idx = 1:num_elements
        element_faces_list = mymesh.element_faces_list(ele_idx,:);
        % go through each face to do (Q,U)'*uhat
        for ll = 1:length(element_faces_list)
            face_idx = element_faces_list(ll);
            temp_uhat = uh_hat ((face_idx-1)*Nuhat+1:face_idx*Nuhat);
            
            Result_matrix(:,ele_idx) = Result_matrix(:,ele_idx)...
                + List_LocSol(:,:,ele_idx,ll)*temp_uhat;
            
        end
        % add (Qwf,Uwf)'
        Result_matrix(:,ele_idx) = Result_matrix(:,ele_idx)...
            +List_LocSol_f(:,ele_idx);
    end
    
    qh = Result_matrix(1:Nq,:);
    uh = Result_matrix(Nq+1:end,:);
    
    
    
end