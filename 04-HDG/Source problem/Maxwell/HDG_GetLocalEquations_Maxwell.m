function [List_LocSol, List_LocSol_f, List_Ns]=HDG_GetLocalEquations_Maxwell(pb, mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,tau_t,tau_n,source_j_1,source_j_2)
    %% ----- Step 0 Set up parameters ---------------------------------------
    
    Nw = (k+1)*(k+2)/2;
    Nu = 2*Nw;
    Np = Nw;
    N_local = Nw+Nu+Np;
    
    Nuhat_t = k+1;
    Nphat = k+1;
    
    N_global = Nuhat_t + Nphat;
    
    num_elements = mymesh.num_elements;
    
    
    %% ----- step 1 get Local Matrices on ref element ----------------------------------
    
    if strcmp(pb,'5')
        % ?time-harmonic Maxwells? equations in
        [Auu,Auur,Auus,Auu3,Buuhat3] = HDG_LocalMatrix(k,GQ1DRef_pts,GQ1DRef_wts);
    else
        error('HDG method for problem type %s has not implemented yet ', pb);
    end
    
    %% ----- step 2 Get Local Solver matrices -------------------------
    
    %%% Local equation
    %
    %          wh                    uhat_t_Fi       0
    %  M_loc* (uh  ) = Sum: N_i * ?          )   +( Proj_j )
    %          ph                    phat_Fi         0
    % >>
    %
    %  wh                  uhat_t_Fi
    %  uh   = Sum: Loc_i *                 + Loc_f
    %  ph                   phat_Fi
    %
    %%%%
    
    List_LocSol_f = zeros(N_local,num_elements,numeric_t); % Local solver for source_j_1 and source_j_2
    
    List_LocSol = zeros(N_local,N_global,num_elements,3,numeric_t); % Local solver for uhat_t and phat
    
    List_Ns = zeros(N_local,N_global,num_elements,3,numeric_t); % helper for Global Eqs
    % go through all the elements to build local solver matrices
    
    for element_idx = 1: num_elements
        
        temp_element = mymesh.element_list(element_idx,:);
        
        vertice_list = mymesh.vertices_list(temp_element(:),:);
        Jk = mymesh.Jacobian_list(element_idx);
        
        uhat_dir_list = mymesh.uhat_dir_list(element_idx,:);
        
        % Local equation matrices
        if strcmp(pb,'5')
            [Ns,M_Loc] = HDG_MaxwellLocalEquations(...  % ^^^^^^^^^^^^^
                                    Jk,vertice_list,tau,...
                                       Auu,Auur,Auus,Auu3,Buuhat3,uhat_dir_list);
        end
        
        % Local solver Q,U
        % [Q;U] = M_loc^-1 * N_i
        
        for jj = 1:3
            List_LocSol(:,:,element_idx,jj) = M_Loc\Ns(:,:,jj);
            % local matrix used for global equations
            % 
            List_Ns(:,:,element_idx,jj) = Ns(:,:,jj);
            List_Ns(Nw+1:Nw+Nu,Nuhat_t+1:end,element_idx,jj)...
                = -1.0 * List_Ns(Nw+1:Nw+Nu,Nuhat_t+1:end,element_idx,jj);
            
        end
            
        % Compute the projection of source_f
        Proj_f = zeros(Nq+Nu,1,numeric_t);
        
        Proj_f(Nq+1:end) = Project_F_to_Wh(Jk,vertice_list,...
            source_f,k,...
            GQ1DRef_pts,GQ1DRef_wts);
        
        % Local solver Qw * f , Uw * f
        % M_loc^-1 * Proj_f
        List_LocSol_f(:,element_idx) = M_Loc\Proj_f;

    end
    
    
end