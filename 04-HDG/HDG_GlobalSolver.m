function [uh,qh,uhat] = HDG_GlobalSolver(pb, mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,tau,source_f,uD,uN)
    
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
        
        % Local equation matrices
        if strcmp(pb,'1')
            [N_1,N_2,N_3,M_Loc] = HDG_PoissionLocalEquations(...
                                    Jk,vertice_list,tau,...
                                       Aqq,Auur,Auus,Auu3,Buuhat3);
        end
        
        % Local solver Q,U
        % M_loc^-1 * N_i
        
        List_LocSol(:,:,1) = M_Loc\N_1;
        List_Ns(:,:,1) = [-N_1(1:Nq,:);N_1(Nq+1:end,:)];
        
        List_LocSol(:,:,2) = M_Loc\N_2;
        List_Ns(:,:,2) = [-N_2(1:Nq,:);N_3(Nq+1:end,:)];
        
        List_LocSol(:,:,3) = M_Loc\N_3;
        List_Ns(:,:,3) = [-N_3(1:Nq,:);N_3(Nq+1:end,:)];
        
        % Compute the projection of source_f
        Proj_f = zeros(Nq+Nu,1,numeric_t);
        
        Proj_f(Nq+1:end) = Project_F_to_Wh(Jk,vertice_list,...
            source_f,k,...
            GQ1DRef_pts,GQ1DRef_wts);
        
        % Local solver Qw * f , Uw * f
        % M_loc^-1 * Proj_f
        List_LocSol_f(:,element_idx) = M_loc\Proj_f;

    end
    
    %% ----- step 3 Assemble Global Matrix and right hand side --------
   
    
    %% ------ Step 4: Solve the Global Equals ------------------------
    
    
    %% ------ Step 5: Recover Local Solutions --------------------------
end