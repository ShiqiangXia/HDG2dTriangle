function test_uh_pts(para)
    
    %% Step 1 : Set up some varibales
   
    pb_type = num2str(para.pb_type);
    Niter = para.Niter;
    
    % bonus parameters
    temp_report_flag = 1;
    plot_log_err_flag = 1;
    posterior_estimate_method = 2;
    latex_table_flag = 1;
    save_flag = 1;
    
    % adaptive 
    mark_flag = 1; % 0: max marking strategy
                   % 1: bulk marking strategy Dorfler , 
                   % 2? equi distribution strategy
                   % 3: fraction marking strategy
                   
    %%%%%% step 1. Set varibales to store results %%%%%%%%%%%%%%%%%%%%%
    err_cal_flag = para.err_cal_flag;
    err_analysis_flag = para.err_analysis_flag;
    mesh_list = zeros(Niter,1,numeric_t);
    tri_list  = zeros(Niter,1,numeric_t); % record # of triangles for each mesh
    
    reduce_ratio = MyParaParse(para.extra_parameters,'reduce_ratio');
    tol_adp = MyParaParse(para.extra_parameters,'tol_adp');
    
    
    if err_cal_flag == 1
        err_uh_list = zeros(Niter,1,numeric_t);
        err_qh_list = zeros(Niter,1,numeric_t);
        
        err_uh_coarse_list = zeros(Niter,1,numeric_t);
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
    
    %%  step 2: Determine what problem we are solving.
    
    if strcmp(pb_type(1),'1')
        %%  Solve PDE problem
        % step 3. Iterative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cprintf('blue','--------------------------------\n')
        cprintf('blue','Start solving PDE problem\n')
        
        coarse_mesh = Build2DMesh(para.structure_flag, para.dom_type,...
                    2*para.h0, ...
                    para.dirichlet_flag, para.neuman_flag, para.geo_parameters{:});
                
        for ii = 1:Niter
            cprintf('blue','Mesh %d ... \n',ii)
            
            % build mesh --------------------------------------------------
            if ii == 1
                mymesh = Build2DMesh(para.structure_flag, para.dom_type,...
                    para.h0, ...
                    para.dirichlet_flag, para.neuman_flag, para.geo_parameters{:});
            else
                if para.refine_flag == 0 % uniform refinement
                   
                    mymesh = mymesh.UniformRefine();
                    
                elseif para.refine_flag == -1 % build a new mesh with h
                    
                    hh = para.h0*(0.5^(ii-1));
                    
                    mymesh = Build2DMesh(para.structure_flag, para.dom_type,...
                        hh, ...
                        para.dirichlet_flag, para.neuman_flag, para.geo_parameters{:});
                    
                else % refine based on marked elements
                    r_f = GetRefineMethod(para.refine_flag); %'R', 'RG', 'NVB'
                    mymesh = mymesh.Refine(marked_elements, r_f);
                end
            end
            
            
            if strcmp(pb_type(2),'0')
                %% Solve source problem
                
                [source_f,uD,uN]=MyParaParse(para.pb_parameters,'source_f','uD','uN');
                [uh,qh,uhat] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau,source_f,uD,uN); 
            end
            
            %% Evalue uh for GQ pts on coarse mesh (merged to squares)
            [GQ_x, GQ_y, hx, hy] = GetPhyGQPts(para.structure_flag,para.dom_type,...
                2*para.h0,GQ1DRef_pts, para.geo_parameters{:});
            
            uh_coarse_GQ_pts = GetUhGQPtsatCoarseMesh(para.order, ...
                coarse_mesh, mymesh, uh, GQ_x, GQ_y);
            
            %%  Calculate function Error ---------------------------------------------
            mesh_list(ii) = GetDof(mymesh, para.order);
            
            if err_cal_flag
                
                [uexact,qexact_1,qexact_2]=MyParaParse(para.pb_parameters,'uexact','qexact_1','qexact_2');
            
                [err_uh_list(ii),err_uh_elewise] = L2Error_scalar(mymesh,uh,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    para.order,uexact);
                
                %PlotElementWiseValue(mymesh,err_uh_elewise,'err-uh elementwise' );

                [err_qh_list(ii),~] = L2Error_vector(mymesh,qh,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    para.order,qexact_1,qexact_2);
                
                err_uh_coarse_list(ii) = L2Error_scalar_Square(uh_coarse_GQ_pts,...
                    uexact, GQ_x, GQ_y, GQ1DRef_wts, hx, hy);
                
                
            end
            
        end
        
        % Step 4. Report reulsts%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if para.report_flag==1 
            ReportProblem(para) 
            % report uh, qh errors
            if  err_cal_flag == 1
                order_uh = GetOrder(mesh_list,err_uh_list);
                order_qh = GetOrder(mesh_list,err_qh_list);
                
                order_uh_coarse = GetOrder(mesh_list,err_uh_coarse_list);
                
                ReportTable('DOF', mesh_list,...
                    'err_uh',err_uh_list,'order', order_uh,...
                    'err_qh',err_qh_list,'order',order_qh,...
                    'err_uh_coar', err_uh_coarse_list, 'order',order_uh_coarse)
            end
        end
        
    end
    

    
end