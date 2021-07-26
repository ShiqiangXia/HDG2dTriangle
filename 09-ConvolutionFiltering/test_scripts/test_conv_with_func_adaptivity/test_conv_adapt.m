function test_conv_adapt(para)
    %----------------------------------------------------------------------
    %% step 1: set up 
    %----------------------------------------------------------------------
   
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
    
    %----------------------------------------------------------------------
    % step 2: Get uH_triangle, vH_triangle and corase mesh TH
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % step 3: Get uh_final, vh_final and fine mesh Th after adaptivity N times
    %----------------------------------------------------------------------
    
    if strcmp(pb_type(1),'1')
        error('Call the wrong script! This script is to test funcitonal adaptivity')
    else
        %%  Solver Functional problem
        pb_text_info = MyParaParse(para.pb_parameters,'pb_text_info');
        
        if strcmp(pb_type(2),'0') % define list to store results
           
            Jh_list = zeros(Niter,1,numeric_t);
            err_Jh_list = zeros(Niter,1,numeric_t);
            Jh_AC_list = zeros(Niter,1,numeric_t);
            err_Jh_AC_list = zeros(Niter,1,numeric_t);
            ACh_list = zeros(Niter,1,numeric_t);
            err_vh_list = zeros(Niter,1,numeric_t);
            estimate_sum_abs_list = zeros(Niter,1,numeric_t);
            if para.post_process_flag == 1
                est_terms_sum_list = zeros(Niter,1,numeric_t);
            end
            
        else
            error('Wrong problem type')
        end
        
        cprintf('blue','--------------------------------\n')
        cprintf('blue','Start solving functional problem\n')

        for ii = 1:Niter
            cprintf('blue','Mesh %d ... \n',ii)
            
            % Build initial mesh
            if ii == 1
                
                mymesh = Build2DMesh(para.structure_flag, para.dom_type,...
                    para.h0, ...
                    para.dirichlet_flag, para.neuman_flag, para.geo_parameters{:});
                
                
                
            else
                if para.refine_flag == 0 % uniform refinement
                    %% uniform refinement
                    
                    mymesh = mymesh.UniformRefine();

                elseif para.refine_flag == -1 % build a new mesh with h
                    %% build a new mesh with a new h
                    
                    hh = para.h0*(0.5^(ii-1));
                    
                    mymesh = Build2DMesh(para.structure_flag, para.dom_type,...
                        hh, ...
                        para.dirichlet_flag, para.neuman_flag, para.geo_parameters{:});
                        
                else % refine based on marked elements
                    %% adpvie refine
                    r_f = GetRefineMethod(para.refine_flag); %'R', 'RG', 'NVB'
                    mymesh = mymesh.Refine(marked_elements, r_f);
 
                end
            end
            
            
            if ii == Niter
                finer_mesh = mymesh;
                finer_mesh.Plot2(0,"Final adaptive mesh")
            elseif ii == 1
                h_coarse = para.h0;
                coarse_mesh = mymesh;
                coarse_mesh.Plot2(0,"Coarse mesh H=" + num2str(h_coarse))
            end
            
            mesh_list(ii) = GetDof(mymesh, para.order);
            tri_list(ii) = mymesh.num_elements;
            if strcmp(pb_type(2),'0') % source problem
                %%  need to solve Primal and Adjoint two problems
                if strcmp(pb_type(3),'1') % Solve Poission source problem
             
                    [source_f,uD,uN]=MyParaParse(para.pb_parameters,'source_f','uD','uN');
                    [source_g,vD,vN]=MyParaParse(para.pb_parameters,'source_g','vD','vN');
                    
                    [uh,qh,uhat] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau,source_f,uD,uN);
                
                    [vh,ph,vhat] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau,source_g,vD,vN);
                    
                    if ii == 1
                        uH_triangle = uh;
                        vH_triangle = vh;
                    end
             
                    if ii == Niter
                        uh_final = uh;
                        vh_final = vh;
                    end
                
                    
                else
                    error('Wrong problem type.')
                end
                
                [Jh,Jh_AC,ACh,ACh_elewise_list,Jh_elewise_list] = LinearFunctional_Elliptic(pb_type(4),pb_type(3),mymesh,...
                                                      uh,qh,uhat,source_f,...
                                                      vh,ph,vhat,source_g,...
                                                      GQ1DRef_pts,GQ1DRef_wts,para.order,para.tau ,0);
                Jh_list(ii) = Jh;
                Jh_AC_list(ii) = Jh_AC;
                ACh_list(ii) = ACh;
                % Calculate error
                if err_cal_flag==1
                    
                    % error of uh, qh
                    [uexact,qexact_1,qexact_2]=MyParaParse(para.pb_parameters,'uexact','qexact_1','qexact_2');

                    [err_uh_list(ii),~] = L2Error_scalar(mymesh,uh,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order,uexact);

                    [err_qh_list(ii),~] = L2Error_vector(mymesh,qh,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order,qexact_1,qexact_2);
                    
                    %PlotElementWiseValue(mymesh,err_uh_elewise,'err-uh elementwise' );

                    [err_Jh_list(ii),err_Jh_AC_list(ii),err_Jh_elewise] ...
                        = Error_Functional(pb_type(4),para.pb_parameters,mymesh,GQ1DRef_pts,GQ1DRef_wts,Jh,Jh_AC,Jh_elewise_list);
                    
                    % error vh, ph;
                    [vexact,pexact_1,pexact_2]=MyParaParse(para.pb_parameters,'vexact','pexact_1','pexact_2');
                    [err_vh_list(ii),err_vh_elewise] = L2Error_scalar(mymesh,vh,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order,vexact);
                    %PlotElementWiseValue(mymesh,err_vh_elewise,'err-vh elementwise' );
                    
                end
                

            end
            % post-processing
            if para.post_process_flag == 1 
 
                if strcmp(pb_type(2),'0') && (strcmp(pb_type(4),'1') ||strcmp(pb_type(4),'2'))
                    %%  poission source problem
                
                    [uhstar,qhstar] = HDG_Local_Postprocess_Elliptic(mymesh,para.order,uh,qh,uhat,para.tau,GQ1DRef_pts,GQ1DRef_wts);
                    [vhstar,phstar] = HDG_Local_Postprocess_Elliptic(mymesh,para.order,vh,ph,vhat,para.tau,GQ1DRef_pts,GQ1DRef_wts);
                    
                    err_uhstar_list(ii) = L2Error_scalar(mymesh,uhstar,...
                        GQ1DRef_pts,GQ1DRef_wts,1,...
                        para.order,uexact);

                    err_qhstar_list(ii)= L2Error_vector(mymesh,qhstar,...
                       GQ1DRef_pts,GQ1DRef_wts,1,...
                        para.order,qexact_1,qexact_2);
                    if posterior_estimate_method == 2
                        [est_terms_sum,est1,est2,est3,est4,est5,pos1,pos2,pos3,pos4]=...
                            Functional_Eh_Estimate_Residual_method(pb_type(4),pb_type(3),mymesh,...
                            uhstar,qhstar,source_f,...
                            vhstar,phstar,source_g,...
                            uh,qh,uhat,...
                            vh,ph,vhat,...
                            GQ1DRef_pts,GQ1DRef_wts,para.order,para.tau,uexact);
                        
                        est_terms_sum_list(ii) = sum(est_terms_sum);
                    end
                end
            end
            
            if para.post_process_flag == 1
                estimate_functinal_elewise = est_terms_sum+ACh_elewise_list;%+est_terms_sum ;%+%;
            else
                estimate_functinal_elewise = ACh_elewise_list;
            end
            
            estimate_sum_abs_list(ii) = sum(abs(estimate_functinal_elewise));
            
            % adaptive mesh refinement
            if para.refine_flag > 0
                percent = MyParaParse(para.extra_parameters,'percent');
                marked_elements = ACh_ErrEstimate(estimate_functinal_elewise,tol_adp,percent,mark_flag);
            end
            
            
        end
 
    %----------------------------------------------------------------------
    % step 4: Get uh_proj, vh_proj by L2 projection of uh, vh to TH
    %         Get uH_square, vH_square by L2 projection of uH_triangle,
    %         vH_triangle to TH
    %----------------------------------------------------------------------
        
        poly_k = para.order;
        N_bd = 2*poly_k;
        spline_degree = poly_k;
        poly_proj = poly_k; 
        y_cut = 5;
        n_level = 3;
        Conv_matrix = Get_convolution_matrix(GQ1DRef_pts,spline_degree,poly_proj+1,GQ1DRef_pts,GQ1DRef_wts);
        
        [GQ_x, GQ_y, hx, hy,Nx_coarse,Ny_coarse] = GetPhyGQPts(para.structure_flag,para.dom_type,...
                h_coarse,GQ1DRef_pts, para.geo_parameters{:});
            
        % project uH_triagnle to uH_square
        uH_square = GetUhL2ProjectionCoarseMesh(poly_k,coarse_mesh,coarse_mesh,...
                uH_triangle,GQ1DRef_pts, GQ1DRef_wts,hx,hy);
        
        % project uh_final to TH
        uh_proj = GetUhL2ProjectionCoarseMesh(poly_k,coarse_mesh,finer_mesh,...
                uh_final,GQ1DRef_pts, GQ1DRef_wts,hx,hy);
            
        uH_square_GQpts = GetUhProjGQpts(uH_square,poly_proj,GQ1DRef_pts);    
        uh_proj_GQpts = GetUhProjGQpts(uh_proj,poly_proj,GQ1DRef_pts);
       
        
        % Convolution Postprocessing
        MH = ConvolutionFiltering(para.dom_type,spline_degree,...
                uH_square, Nx_coarse, Ny_coarse, N_bd, Conv_matrix);
            
        Mh_proj = ConvolutionFiltering(para.dom_type,spline_degree,...
                uh_proj, Nx_coarse, Ny_coarse, N_bd, Conv_matrix);
        % some plots
        flag_2D_plot = 0;
        if flag_2D_plot == 1
            
            
            Plot2D(para.dom_type, GQ_x, GQ_y, MH, "$u_H^*$ on coarse mesh ",0,"")
           
            Plot2D(para.dom_type, GQ_x, GQ_y, Mh_proj, "$u_{h,proj}^*$ on coarse mesh",0,"")

            
        end
            
        if err_cal_flag
            
            uexact_GQ_pts = GetUexactGQpts(uexact, GQ_x, GQ_y);
            Plot2D(para.dom_type, GQ_x, GQ_y, uexact_GQ_pts-uH_square_GQpts, "$u_{H,square}$ on coarse mesh ",0,"")
            Plot2D(para.dom_type, GQ_x, GQ_y, uexact_GQ_pts-uh_proj_GQpts, "$u_{h,proj}$ on coarse mesh ",0,"")
            
            
            uexact_GQ_pts = RemoveBdElements(para.dom_type, uexact_GQ_pts,N_bd,Nx_coarse, Ny_coarse);
            
            err_uHstar = L2Error_scalar_Square(MH,...
                    uexact_GQ_pts, GQ1DRef_wts, hx, hy);
                
            err_uhstar = L2Error_scalar_Square(Mh_proj,...
                    uexact_GQ_pts, GQ1DRef_wts, hx, hy);
                
            flag_plot_diff = 1;
            if flag_plot_diff == 1
                save_flag = 0;
                name_text = "nothing";
                
                Plot2D(para.dom_type, GQ_x, GQ_y, uexact_GQ_pts-MH,...
                    "$u- u_H^*$",save_flag,...
                    name_text)

                Plot2D(para.dom_type, GQ_x, GQ_y, uexact_GQ_pts-Mh_proj,...
                    "$u - u_{h,proj}^*$ ",save_flag,...
                    name_text) 

            end
                
        end
        
        fprintf('err_uH* : %.2e\n',err_uHstar)
        fprintf('err_uh*_proj: %.2e\n',err_uhstar);
        
        
    end
    
   
    
   
    
    %----------------------------------------------------------------------
    % step 5: Do convolution and compare uH* and uH_proj*
    %----------------------------------------------------------------------
end