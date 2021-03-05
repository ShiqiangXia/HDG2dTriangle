function ProblemDriver(para)
    % This is the main problem driver.
    pb_type = num2str(para.pb_type);
    Niter = para.Niter;
    
    %%%%%% step 1. Set varibales to store results %%%%%%%%%%%%%%%%%%%%%
        
    mesh_list = zeros(Niter,1,numeric_t);
    err_uh_list = zeros(Niter,1,numeric_t);
    err_qh_list = zeros(Niter,1,numeric_t);
    
    if strcmp(pb_type(2),'1')
        % eigenvalue problem
        % maybe use matlab inputParser later
        [Neig,Max_iter,Tol_eig] = MyParaParse(para.extra_parameters,...
            'Neig','Max_iter','tol_eig');
    %             extra_para = para.extra_parameters;
    %             extra_para = reshape(extra_para,[],2)';
    %             Neig = extra_para{find(strcmp(extra_para,'Neig')),2};
        err_lamh_list = zeros(Niter,Neig,numeric_t);
        lamh_list = zeros(Niter,Neig,numeric_t);
    end
    
    [GQ1DRef_pts,GQ1DRef_wts] = GaussQuad(para.GQ_deg);
    
    
    if para.post_process_flag == 1
        err_uhstar_list = zeros(Niter,1,numeric_t);
        err_qhstar_list = zeros(Niter,1,numeric_t);
    end
    
    % step 2: Determine what problem we are solving.
    
    if strcmp(pb_type(1),'1')
        
        % Solve PDE problem
    
    
        %%%%%%% step 3. Iterative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ii = 1:Niter
            
            % build mesh --------------------------------------------------
            if ii == 1
                mymesh = Build2DMesh(para.structure_flag, para.dom_type,...
                    para.h0, ...
                    para.dirichlet_flag, para.neuman_flag, para.geo_parameters);
            else
                if para.refine_flag == 0 % uniform refinement
                    %hh = para.h0*(0.5^(ii-1));
                    %
                    mymesh = mymesh.UniformRefine();
                    
                else % refine based on marked elements
                    r_f = 'RGB'; %'R', 'RG', 'NVB'
                    mymesh = mymesh.Refine(marked_elements, r_f);
                end
            end
            
            mesh_list(ii) = GetDof(mymesh, para.order);
            
            % -------------------------------------------------------------
                       
            % Solve -------------------------------------------------------
            
            
            if strcmp(pb_type(3),'1') && strcmp(pb_type(2),'0')
                % Solve Poission source problem
                
                [source_f,uD,uN]=MyParaParse(para.pb_parameters,'source_f','uD','uN');
                [uh,qh,uhat] = HDG_Poission(mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau,source_f,uD,uN);
                
            % -------------------------------------------------------------
            
            
            elseif  strcmp(pb_type(3),'1') && strcmp(pb_type(2),'1')
                % Solve Poission eigen problem
                
                [lamh,uh,qh,uhat] = HDG_PoissionEig(mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                    para.order,para.tau, Neig,Max_iter,Tol_eig);
                lamh_list(ii) = lamh;
                err_lamh_list(ii) = EigenError(para,lamh);
            % -------------------------------------------------------------
            else
                error('pb type not implemented yet')
            end 
            
            % -------------------------------------------------------------
            
            
            % Posterior error estimate if needed---------------------------
            if para.refine_flag == 1 && strcmp(pb_type(2),'0')
                [Tol_adp,Percent] = MyParaParse(para.extra_parameters,'Tol_adp','Percent');
                marked_elements = ErrorEstimate(mymesh,uh,qh,uhat,source_f,...
                    Tol_adp,Percent);
            end
            % -------------------------------------------------------------
            
            % Calculate Error ---------------------------------------------
            [uexact,qexact]=MyParaParse(para.pb_parameters,'uexact','qexact');
            
            err_uh_list(ii) = L2Error_scalar(mymesh,uh,...
                GQ1DRef_pts,GQ1DRef_wts,0,...
                para.order,uexact);

            err_qh_list(ii) = L2Error_vector(mymesh,qh,...
                GQ1DRef_pts,GQ1DRef_wts,0,...
                para.order,qexact);
            
            % post processed error 
            if para.post_process_flag == 1
                [uhstar,quhstar] = HDG_Local_Postprocess(mymesh,para.order,uh,qh,uhat);

                err_uhstar_list(ii) = L2Error_scalar(mymesh,uhstar,...
                    GQ1DRef_pts,GQ1DRef_wts,1,...
                para.order,para.pb_parameters);

                err_qhstar_list(ii)= L2Error_vector(mymesh,quhstar,...
                   GQ1DRef_pts,GQ1DRef_wts,1,...
                para.order,para.pb_parameters);
            end
            % ------------------------------------------------------------- 
            
            
            % visualization -----------------------------------------------
            if ii == Niter && para.visualize_flag==1
                basis_flag = 0;
                Plot(mymesh,uh,para.order, GQ1DRef_pts,basis_flag );
            end
            % ------------------------------------------------------------- 
            
            
        end
        
        % Step 4. Report reulsts%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if para.report_flag==1
            
            ReportProblem(para) 
            
            ReportTable('dof', mesh_list,...
                'err_uh',err_uh_list,...
                'err_qh',err_qh_list)
            
            if para.post_process_flag == 1
                ReportTable('dof', mesh_list,...
                    'err_uhstar',err_uhstar_list,...
                    'err_qhstar',err_qhstar_list)
            end
            
            if strcmp(pb_type(2),'1')
                ReportTable('dof', mesh_list,...
                    'lamh',err_lamh_list)
            end
            
        end
        % -----------------------------------------------------------------
    
    elseif strcmp(pb_type(1),'2')
        
        % Solver Functional problem
        if strcmp(pb_type(2),'0') % source problem
            Jh_list = zeros(Niter,1,numeric_t);
            err_Jh_list = zeros(Niter,1,numeric_t);
            Jh_AC_list = zeros(Niter,1,numeric_t);
            err_Jh_AC_list = zeros(Niter,1,numeric_t);
            ACh_list = zeros(Niter,1,numeric_t);
            
        elseif strcmp(pb_type(2),'1') % eigenproblem
            err_lamh_AC_list = zeros(Niter,1,numeric_t);
            lamh_AC_list = zeros(Niter,1,numeric_t);
            lamh2_list = zeros(Niter,1,numeric_t);
            ACh_list = zeros(Niter,1,numeric_t);
        else
            error('Wrong problem type')
        end
        
        for ii = 1:Niter
            
            % build mesh --------------------------------------------------
            if ii == 1
                mymesh = Build2DMesh(para.structure_flag, para.dom_type,...
                    para.h0, ...
                    para.dirichlet_flag, para.neuman_flag, para.geo_parameters);
            else
                if para.refine_flag == 0 % uniform refinement
                    %hh = para.h0*(0.5^(ii-1));
                    %
                    mymesh = mymesh.UniformRefine();
                    
                else % refine based on marked elements
                    r_f = 'RGB'; %'R', 'RG', 'NVB'
                    mymesh = mymesh.Refine(marked_elements, r_f);
                end
            end
            % -------------------------------------------------------------
            mesh_list(ii) = GetDof(mymesh, para.order);
            
            % Solve -------------------------------------------------------
            if strcmp(pb_type(2),'0') % source problem
                % need to solve Primal and Adjoint two problems
                if strcmp(pb_type(3),'1') % Poission source problem
             
                    [source_f,uD,uN]=MyParaParse(para.pb_parameters,'source_f','uD','uN');
                    [source_g,vD,vN]=MyParaParse(para.pb_parameters,'source_g','vD','vN');
                    
                    [uh,qh,uhat] = HDG_Poission(mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau,source_f,uD,uN);
                
                    [vh,ph,vhat] = HDG_Poission(mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau,source_g,vD,vN);
                    
                else
                    error('Wrong problem type.')
                end
                
                [Jh,Jh_AC,ACh,ACh_elewise_list] = LinearFunctional(pb_type(4),mymesh,...
                                                      uh,qh,uhat,source_f,uD,uN,...
                                                      vh,ph,vhat,source_g,vD,vN);
                Jh_list(ii) = Jh;
                Jh_AC_list(ii) = Jh_AC;
                ACh_list(ii) = ACh;

                % CalError ------------------------------------------------
                [uexact,qexact]=MyParaParse(para.pb_parameters,'uexact','qexact');
            
                err_uh_list(ii) = L2Error_scalar(mymesh,uh,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    para.order,uexact);

                err_qh_list(ii) = L2Error_vector(mymesh,qh,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    para.order,qexact);
                
                [err_Jh_list(ii),err_Jh_AC_list(ii)] ...
                    = Error_Functional(Jh,Jh_AC,para);

            elseif strcmp(pb_type(2),'1') % eigenproblem
                % only need to solve one eigenvlaue problem
                if strcmp(pb_type(3),'1') % poission eigen problem
                    
                    [lamh,uh,qh,uhat] = HDG_PoissionEig(mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                    para.order,para.tau, Neig,Max_iter,Tol_eig);
                    lamh_list(ii) = lamh;
                    err_lamh_list(ii) = EigenError(para,lamh);
                    
                end
                
                [lamh2, lamh_AC,ACh,ACh_elewise_list]=EigenvalueFunctional(mymesh,uh,qh,uhat); 
                lamh2_list(ii) = lamh2;
                lamh_AC_list(ii) = lamh_AC;
                ACh_list(ii) = ACh;
                
                [err_lamh_list(ii),err_lamh_AC_list(ii)] ...
                    = Error_Functional(lamh2,lamh_AC,para);
                 
                
            end
            
            % Posterior error estimate if needed--------------------------- 
            if para.refine_flag == 1
                [Tol_adp,Percent] = MyParaParse(para.extra_parameters,'Tol_adp','Percent');
                marked_elements = ACh_ErrEstimate(ACh_elewise_list,Tol_adp,Percent);
            end
            % -------------------------------------------------------------
               
        end
        
        % Report reulsts---------------------------------------------------
        if para.report_flag==1
            
            ReportProblem(para) 
            
            if strcmp(pb_type(2),'0')
                ReportTable('dof', mesh_list,...
                    'err_uh',err_uh_list,...
                    'err_qh',err_qh_list)
            
                ReportTable('dof', mesh_list,...
                    'err_Jh',err_Jh_list,...
                    'ACh',ACh_list,...
                    'err_Jh_AC',err_Jh_AC_list)
                
            elseif strcmp(pb_type(2),'1')
                ReportTable('dof', mesh_list,...
                    'err_lamh',err_lamh_list,...
                    'ACh',ACh_list,...
                    'err_lamh_AC',err_lamh_AC_list)
                
            else
                
            end
            
            
        end
               
        
    else
        error('Pb type is not incorrect, please double check and see Parameter obj')
    end
    
end