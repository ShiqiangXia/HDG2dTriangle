function ProblemDriver(para)
    % This is the main problem driver.
    pb_type = num2str(para.pb_type);
    Niter = para.Niter;
    
    %%%%%% step 1. Set varibales to store results %%%%%%%%%%%%%%%%%%%%%
    err_cal_flag = para.err_cal_flag;
    if err_cal_flag == 1
        mesh_list = zeros(Niter,1,numeric_t);
        err_uh_list = zeros(Niter,1,numeric_t);
        err_qh_list = zeros(Niter,1,numeric_t);
    end
    
    if strcmp(pb_type(2),'1')
        % eigenvalue problem
        % maybe use matlab inputParser later
        [Neig,Max_iter,Tol_eig] = MyParaParse(para.extra_parameters,...
            'Neig','Max_iter','tol_eig');
        lamh_list = zeros(Niter,Neig,numeric_t);
        
        err_lamh_list = zeros(Niter,Neig,numeric_t);
        
    end
    
    [GQ1DRef_pts,GQ1DRef_wts] = GaussQuad(para.GQ_deg);
    
    
    if para.post_process_flag == 1 && err_cal_flag == 1
        err_uhstar_list = zeros(Niter,1,numeric_t);
        err_qhstar_list = zeros(Niter,1,numeric_t);
    end
    
    % step 2: Determine what problem we are solving.
    
    if strcmp(pb_type(1),'1')
        
        % Solve PDE problem
    
    
        %%%%%%% step 3. Iterative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cprintf('blue','--------------------------------\n')
        cprintf('blue','Starting solving PDE problem\n')
        for ii = 1:Niter
            cprintf('blue','Mesh %d ... \n',ii)
            % build mesh --------------------------------------------------
            if ii == 1
                mymesh = Build2DMesh(para.structure_flag, para.dom_type,...
                    para.h0, ...
                    para.dirichlet_flag, para.neuman_flag, para.geo_parameters{:});
            else
                if para.refine_flag == 0 % uniform refinement
                    hh = para.h0*(0.5^(ii-1));
                    %
                    mymesh = Build2DMesh(para.structure_flag, para.dom_type,...
                        hh, ...
                        para.dirichlet_flag, para.neuman_flag, para.geo_parameters{:});
                
                    %mymesh = mymesh.UniformRefine();
                    
                else % refine based on marked elements
                    r_f = GetRefineMethod(para.refine_flag); %'R', 'RG', 'NVB'
                    mymesh = mymesh.Refine(marked_elements, r_f);
                end
            end
            
            
            
            % -------------------------------------------------------------
                       
            % Solve -------------------------------------------------------
            
            
            if strcmp(pb_type(3),'1') && strcmp(pb_type(2),'0')
                % Solve Poission source problem
                
                [source_f,uD,uN]=MyParaParse(para.pb_parameters,'source_f','uD','uN');
                [uh,qh,uhat] = HDG_GlobalSolver(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau,source_f,uD,uN);
                
            % -------------------------------------------------------------
            
            
            elseif  strcmp(pb_type(3),'1') && strcmp(pb_type(2),'1')
                % Solve Poission eigen problem
                [lamh,uh_Neig,qh_Neig,uhat_Neig] = HDG_PoissionEig(mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                    para.order,para.tau, Neig,Max_iter,Tol_eig);
                lamh_list(ii) = lamh;
                
                
                err_lamh_list(ii) = EigenError(para,lamh);
          
                if err_cal_flag
                    [tag_eig] = MyParaParse(para.pb_parameters,'tag_eig');
                    uh = uh_Neig(:,:,:,tag_eig);
                    qh = qh_Neig(:,:,:,tag_eig);
                    uhat = uhat_Neig(:,:,:,tag_eig);
                end
            % -------------------------------------------------------------
            else
                error('pb type not implemented yet')
            end 
            
            % -------------------------------------------------------------
            
            
            % Posterior error estimate if needed---------------------------
            if para.refine_flag == 1 && strcmp(pb_type(2),'0')
                [tol_adp,percent] = MyParaParse(para.extra_parameters,'tol_adp','percent');
                marked_elements = ErrorEstimate(mymesh,uh,qh,uhat,source_f,...
                    tol_adp,percent);
            end
            % -------------------------------------------------------------
            
            % Calculate Error ---------------------------------------------
            
            
            if err_cal_flag
                
                mesh_list(ii) = GetDof(mymesh, para.order);
                
                [uexact,qexact_1,qexact_2]=MyParaParse(para.pb_parameters,'uexact','qexact_1','qexact_2');
            
                [err_uh_list(ii),~] = L2Error_scalar(mymesh,uh,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    para.order,uexact);

                [err_qh_list(ii),~] = L2Error_vector(mymesh,qh,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    para.order,qexact_1,qexact_2);
            end
            
            % post processed error 
            if para.post_process_flag == 1
                [uhstar,quhstar] = HDG_Local_Postprocess(mymesh,para.order,uh,qh,uhat);
                
                if err_cal_flag == 1
                    err_uhstar_list(ii) = L2Error_scalar(mymesh,uhstar,...
                        GQ1DRef_pts,GQ1DRef_wts,1,...
                    para.order,para.pb_parameters);

                    err_qhstar_list(ii)= L2Error_vector(mymesh,quhstar,...
                       GQ1DRef_pts,GQ1DRef_wts,1,...
                    para.order,para.pb_parameters);
                end
            end
            % ------------------------------------------------------------- 
            
            
            % visualization -----------------------------------------------
            if ii == Niter && para.visualize_flag==1
                basis_flag = 0;
                My2DTriPlot(mymesh,uh,para.order, GQ1DRef_pts,basis_flag );
            end
            % ------------------------------------------------------------- 
            
            
        end
        
        % Step 4. Report reulsts%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if para.report_flag==1 
            
            ReportProblem(para) 
            
            if strcmp(pb_type(2),'1')
                ReportTable('dof', mesh_list,...
                    'lamh',err_lamh_list)
            end
            
            if  err_cal_flag == 1
                ReportTable('dof', mesh_list,...
                    'err_uh',err_uh_list,...
                    'err_qh',err_qh_list)

                if para.post_process_flag == 1
                    ReportTable('dof', mesh_list,...
                        'err_uhstar',err_uhstar_list,...
                        'err_qhstar',err_qhstar_list)
                end
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
                    para.dirichlet_flag, para.neuman_flag, para.geo_parameters{:});
            else
                if para.refine_flag == 0 % uniform refinement
                    %hh = para.h0*(0.5^(ii-1));
                    %
                    mymesh = mymesh.UniformRefine();
                    
                else % refine based on marked elements
                    r_f = GetRefineMethod(para.refine_flag); %'R', 'RG', 'NVB'
                    mymesh = mymesh.Refine(marked_elements, r_f);
                end
            end
            % -------------------------------------------------------------
            
            
            % Solve -------------------------------------------------------
            if strcmp(pb_type(2),'0') % source problem
                % need to solve Primal and Adjoint two problems
                if strcmp(pb_type(3),'1') % Poission source problem
             
                    [source_f,uD,uN]=MyParaParse(para.pb_parameters,'source_f','uD','uN');
                    [source_g,vD,vN]=MyParaParse(para.pb_parameters,'source_g','vD','vN');
                    
                    [uh,qh,uhat] = HDG_GlobalSolver(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau,source_f,uD,uN);
                
                    [vh,ph,vhat] = HDG_GlobalSolver(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau,source_g,vD,vN);
                    
                else
                    error('Wrong problem type.')
                end
                
                [Jh,Jh_AC,ACh,ACh_elewise_list] = LinearFunctional(pb_type(4),mymesh,...
                                                      uh,qh,uhat,source_f,uD,uN,...
                                                      vh,ph,vhat,source_g,vD,vN,...
                                                      GQ1DRef_pts,GQ1DRef_wts,para.order);
                Jh_list(ii) = Jh;
                Jh_AC_list(ii) = Jh_AC;
                ACh_list(ii) = ACh;

                % CalError ------------------------------------------------
                
                
                if err_cal_flag==1
                    mesh_list(ii) = GetDof(mymesh, para.order);
                    
                    [uexact,qexact_1,qexact_2]=MyParaParse(para.pb_parameters,'uexact','qexact_1','qexact_2');

                    [err_uh_list(ii),~] = L2Error_scalar(mymesh,uh,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order,uexact);

                    [err_qh_list(ii),~] = L2Error_vector(mymesh,qh,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order,qexact_1,qexact_2);

                    [err_Jh_list(ii),err_Jh_AC_list(ii)] ...
                        = Error_Functional(Jh,Jh_AC,para);
                end

            elseif strcmp(pb_type(2),'1') % eigenproblem
                % only need to solve one eigenvlaue problem
                if strcmp(pb_type(3),'1') % poission eigen problem
                    
                    [lamh,uh_Neig,qh_Neig,uhat_Neig] = HDG_PoissionEig(mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                    para.order,para.tau, Neig,Max_iter,Tol_eig);
                    lamh_list(ii) = lamh;
                    err_lamh_list(ii) = EigenError(para,lamh);
                    
                end
                
                [lamh2,lamh_AC,ACh,ACh_elewise_list]=...
                    EigenvalueFunctional(mymesh,...
                    uh_Neig,qh_Neig,uhat_Neig,...
                    GQ1DRef_pts,GQ1DRef_wts,para.order); 
                
                lamh2_list(ii) = lamh2;
                lamh_AC_list(ii) = lamh_AC;
                ACh_list(ii) = ACh;
                
                if err_cal_flag == 1
                    
                    [err_lamh_list(ii),err_lamh_AC_list(ii)] ...
                        = Error_Functional(lamh2,lamh_AC,para);
                end
                 
                
            end
            
            % Posterior error estimate if needed--------------------------- 
            if para.refine_flag == 1
                [tol_adp,percent] = MyParaParse(para.extra_parameters,'tol_adp','percent');
                marked_elements = ACh_ErrEstimate(ACh_elewise_list,tol_adp,percent);
            end
            % -------------------------------------------------------------
               
        end
        
        % Report reulsts---------------------------------------------------
        if para.report_flag==1 
            
            ReportProblem(para) 
            
            if strcmp(pb_type(2),'0')
                
                ReportTable('dof', mesh_list,...
                    'err_Jh',err_Jh_list,...
                    'ACh',ACh_list,...
                    'err_Jh_AC',err_Jh_AC_list)
                
                if err_cal_flag==1
                    ReportTable('dof', mesh_list,...
                        'err_uh',err_uh_list,...
                        'err_qh',err_qh_list)
                end
           
                
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