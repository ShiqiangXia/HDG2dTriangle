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
        err_uh_proj_coarse_list = zeros(Niter,1,numeric_t);
        err_uhstar_coarse_list = zeros(Niter,1,numeric_t);
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
        
        h_corase = 2*para.h0;
        coarse_mesh = Build2DMesh(para.structure_flag, para.dom_type,...
                    h_corase, ...
                    para.dirichlet_flag, para.neuman_flag, para.geo_parameters{:});
                
        N_bd = 0;%2*poly_k;
        
        poly_k = para.order;
        spline_degree = poly_k;
        poly_proj = poly_k;  
        
        %LGL_pts = JacobiGL(0,0,2*poly_k+1);% 2k+2 Gauss Lobato points for polynomial of degree 2k+1
        Conv_matrix = Get_convolution_matrix(GQ1DRef_pts,spline_degree,poly_proj+1,GQ1DRef_pts,GQ1DRef_wts);
        
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
%             h_corase = para.h0*(0.5^(ii-1));
%             coarse_mesh = mymesh;
%             
            
            if strcmp(pb_type(2),'0')
                %% Solve source problem
                
                [source_f,uD,uN]=MyParaParse(para.pb_parameters,'source_f','uD','uN');
                [uh,qh,uhat] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    poly_k, para.tau,source_f,uD,uN); 
            end
            
            %% Evalue uh for GQ pts on coarse mesh (merged to squares)
            [GQ_x, GQ_y, hx, hy,Nx_coarse,Ny_coarse] = GetPhyGQPts(para.structure_flag,para.dom_type,...
                h_corase,GQ1DRef_pts, para.geo_parameters{:});
            
            uh_coarse_GQ_pts = GetUhGQPtsatCoarseMesh(poly_k, ...
                coarse_mesh, mymesh, uh, GQ_x, GQ_y);
            PlotUhcut(uh_coarse_GQ_pts, hx,10,3, GQ_x,'uh-coarse','uh-corase on coarse mesh' )
            
            
            %%%%%%%%%% project uh to P_k on coarse mesh and then postprocessing
            
            uh_coarse_proj = GetUhProjCoarseMesh(poly_proj,uh_coarse_GQ_pts,GQ1DRef_pts);
            
            uh_proj_GQpts = GetUhProjGQpts(uh_coarse_proj,poly_proj,GQ1DRef_pts);
            PlotUhcut(uh_proj_GQpts, hx,10,3, GQ_x,'uh-proj','uh proj on coarse mesh' )
            

            M = ConvolutionFiltering(para.dom_type,spline_degree,...
                uh_coarse_proj, Nx_coarse, Ny_coarse, N_bd, Conv_matrix);
            
            PlotUhcut(M, hx,10,3, GQ_x,'uh*','uh* on coarse mesh' )
            %M = ConvertUhPts(M,poly_k,LGL_pts,GQ1DRef_pts);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            

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
                
                uexact_GQ_pts = GetUexactGQpts(uexact, GQ_x, GQ_y);
                
                err_uh_coarse_list(ii) = L2Error_scalar_Square(uh_coarse_GQ_pts,...
                    uexact_GQ_pts, GQ1DRef_wts, hx, hy);
                
                err_uh_proj_coarse_list(ii)=L2Error_scalar_Square(uh_proj_GQpts,...
                    uexact_GQ_pts, GQ1DRef_wts, hx, hy);
                
                
                err_uhstar_coarse_list(ii) = L2Error_scalar_Square(M,...
                    uexact_GQ_pts, GQ1DRef_wts, hx, hy);
                
                % Plot
%                 PlotUhcut(uexact_GQ_pts-uh_coarse_GQ_pts, hx,10,3, GQ_x,'u - uh-coarse','u - uh-coarse on coarse mesh' )
%                 PlotUhcut(uexact_GQ_pts-uh_proj_GQpts, hx,10,3, GQ_x,'u - uh-proj','u - uh-proj on coarse mesh' )
%                 PlotUhcut(uexact_GQ_pts-M, hx,10,3, GQ_x,'u - uh*','u - uh* on coarse mesh' )
%                 
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
                order_uh_proj_coarse = GetOrder(mesh_list,err_uh_proj_coarse_list);
                order_uhstar_coarse = GetOrder(mesh_list,err_uhstar_coarse_list);
                
                ReportTable('DOF', mesh_list,...
                    'err_uh',err_uh_list,'order', order_uh,...
                    'err_qh',err_qh_list,'order',order_qh)
                 ReportTable('DOF', mesh_list,...
                    'err_uh_coar', err_uh_coarse_list, 'order',order_uh_coarse,...
                    'err_uh_proj_coar', err_uh_proj_coarse_list, 'order',order_uh_proj_coarse,...
                    'err_uhstar', err_uhstar_coarse_list, 'order',order_uhstar_coarse)
            end
        end
        
    end
    

    
end