function test_conv_adapt_with_extraction(para,Ncoarse, N_outer_adap_steps )
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
    save_flag = 0;
    
    % adaptive 
    mark_flag = 1; % 0: max marking strategy
                   % 1: bulk marking strategy Dorfler , 
                   % 2? equi distribution strategy
                   % 3: fraction marking strategy
                   
    %%%%%% step 1. Set varibales to store results %%%%%%%%%%%%%%%%%%%%%
    
    flag_conv = 1;
    
    err_cal_flag = para.err_cal_flag;
    err_analysis_flag = para.err_analysis_flag;
    mesh_list = zeros(Niter,1,numeric_t);

    tri_list  = zeros(Niter,1,numeric_t); % record # of triangles for each mesh
    h0_list = zeros(Ncoarse,1,numeric_t);
    mesh_comp_list = zeros(Ncoarse,1,numeric_t);
    
    reduce_ratio = MyParaParse(para.extra_parameters,'reduce_ratio');
    tol_adp = MyParaParse(para.extra_parameters,'tol_adp');
    
    
    if err_cal_flag == 1
        err_uh_list = zeros(Niter,1,numeric_t);
        err_qh_list = zeros(Niter,1,numeric_t);
        err_uH_list = zeros(Ncoarse,1,numeric_t);
        err_uH_star_list = zeros(Ncoarse,1,numeric_t);
        err_uh_proj_list = zeros(Ncoarse,1,numeric_t);
        err_uh_proj_star_list = zeros(Ncoarse,1,numeric_t);
        err_uH_star_comp_list = zeros(Ncoarse,1,numeric_t);
        err_Jh_coarse_list = zeros(Ncoarse,1,numeric_t);
        err_Jh_AC_coarse_list = zeros(Ncoarse,1,numeric_t);
        ACh_coarse_list = zeros(Ncoarse,1,numeric_t);
        
        err_Jh_star_list = zeros(Ncoarse,1,numeric_t);
        ACh_star_list = zeros(Ncoarse,1,numeric_t);
        err_Jh_ACh_star_list = zeros(Ncoarse,1,numeric_t);
        
    end
    
    if strcmp(pb_type(2),'0') 
        % define list to store results
        Jh_list = zeros(Niter,1,numeric_t);
        err_Jh_fine_list = zeros(Niter,1,numeric_t);
        Jh_AC_list = zeros(Niter,1,numeric_t);
        err_Jh_AC_fine_list = zeros(Niter,1,numeric_t);
        ACh_list = zeros(Niter,1,numeric_t);
        err_vh_list = zeros(Niter,1,numeric_t);
        estimate_sum_abs_list = zeros(Niter,1,numeric_t);
        if para.post_process_flag == 1
            est_terms_sum_list = zeros(Niter,1,numeric_t);
        end
        if strcmp(pb_type(3),'1')
            % get problem data
            [source_f,uD,uN]=MyParaParse(para.pb_parameters,'source_f','uD','uN');
            [source_g,vD,vN]=MyParaParse(para.pb_parameters,'source_g','vD','vN');
            [uexact,qexact_1,qexact_2]=MyParaParse(para.pb_parameters,'uexact','qexact_1','qexact_2');
            [vexact,pexact_1,pexact_2]=MyParaParse(para.pb_parameters,'vexact','pexact_1','pexact_2');
        end
    elseif strcmp(pb_type(2),'1')
        % eigenvalue problem
        
        [Neig,Max_iter,Tol_eig] = MyParaParse(para.extra_parameters,...
            'Neig','Max_iter','tol_eig');
        lamh_list = zeros(Niter,Neig,numeric_t);
        err_lamh_list = zeros(Niter,Neig,numeric_t);
        
    else
        error('Wrong problem type')
    end
    
  
    
    [GQ1DRef_pts,GQ1DRef_wts] = GaussQuad(para.GQ_deg);
    
    
    if para.post_process_flag == 1 && err_cal_flag == 1
        err_uhstar_list = zeros(Niter,1,numeric_t);
        err_qhstar_list = zeros(Niter,1,numeric_t);

    end
    % precompute the convolution matrix
    if flag_conv == 1
        poly_k = para.order;  
        poly_outer_k = 2*poly_k;
        spline_degree = poly_k;
        poly_proj = poly_k;
        
        y_cut = 5;
        n_level = 3;
        
        %order1_dist = 0.3;
        
        if poly_k ==1
            N_bd =2*poly_k;%0.2/hx; % 0.05
        elseif poly_k == 2
            N_bd =2*poly_k;%0.2/hx;
        else
            N_bd =2*poly_k; %0.3/hx;
        end
        % precompute the convolution matrix
        
        Conv_matrix = Get_convolution_matrix(GQ1DRef_pts,spline_degree,...
            poly_proj+1,GQ1DRef_pts,GQ1DRef_wts);
        
        Conv_matrix_local_post = Get_convolution_matrix(GQ1DRef_pts,spline_degree,...
            poly_proj+2,GQ1DRef_pts,GQ1DRef_wts);
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
        
        
        
        cprintf('blue','--------------------------------\n')
        cprintf('blue','Start solving functional problem\n')
        
        h_coarse_0 =para.h0;
        
        for jj = 1: Ncoarse
            % uniform refine of background coarse mesh
            
            h0 = h_coarse_0*(0.5^(jj-1));
            %h0 = h_coarse_0 -(jj-1)*0.02;
            h0_list(jj) = h0;
            cprintf('blue','\n--------------------------------\n')
            cprintf('blue','Coarse mesh %d, H=%.2e \n\n',jj,h0);
            cprintf('blue','adaptive steps \n');
            
            for ii = 1:Niter
                
                cprintf('blue','Mesh %d ... \n',ii)

                % Build initial mesh
                if ii == 1
                    mymesh = Build2DMesh(para.structure_flag, para.dom_type,...
                        h0, ...
                        para.dirichlet_flag, para.neuman_flag, para.geo_parameters{:});
                else
                    if para.refine_flag == 0 % uniform refinement
                        %% uniform refinement

                        mymesh = mymesh.UniformRefine();

                    elseif para.refine_flag == -1 % build a new mesh with h
                        %% build a new mesh with a new h

                        hh = h0*(0.5^(ii-1));

                        mymesh = Build2DMesh(para.structure_flag, para.dom_type,...
                            hh, ...
                            para.dirichlet_flag, para.neuman_flag, para.geo_parameters{:});

                    else % refine based on marked elements
                        %% adpvie refine
                        r_f = GetRefineMethod(para.refine_flag); %'R', 'RG', 'NVB'
                        mymesh = mymesh.Refine(marked_elements, r_f);

                    end
                end

                flag_mesh_plot = 0;
                if ii == Niter
                    finer_mesh = mymesh;
                    if flag_mesh_plot==1
                        finer_mesh.Plot2(0,"Final adaptive mesh" + num2str(ii) + " based on coarse mesh H = " + num2str(h0))
                        file_name = "k"+num2str(para.order)+"_Adapt_Mesh"+ num2str(ii);
                        if save_flag == 1
                            savefig(gcf,file_name);
                        end
                    end
                end
                if ii == 1
                    h_coarse = h0;
                    coarse_mesh = mymesh;
                    if flag_mesh_plot==1
                        coarse_mesh.Plot2(0,"Coarse mesh H=" + num2str(h_coarse));
                        file_name = "k"+num2str(para.order)+"_Coarse_Mesh";
                        if save_flag == 1
                            savefig(gcf,file_name);
                        end
                    end
                end

                mesh_list(ii) = GetDof(mymesh, para.order);
                tri_list(ii) = mymesh.num_elements;
                
                if strcmp(pb_type(2),'0') % source problem
                    %%  need to solve Primal and Adjoint two problems
                    if strcmp(pb_type(3),'1') % Solve Poission source problem
                        
%                         [uh,qh,uhat] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
%                         para.order, para.tau,source_f,uD,uN);
% 
%                         [vh,ph,vhat] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
%                         para.order, para.tau,source_g,vD,vN);
% 
%                         if ii == 1
%                             uH_triangle = uh;
%                             vH_triangle = vh;
%                         end

%                         if ii == Niter
%                             uh_final = uh;
%                             vh_final = vh;
%                         end


%                     else
%                         error('Wrong problem type.')
%                     end
% 
%                     [Jh,Jh_AC,ACh,ACh_elewise_list,Jh_elewise_list] = LinearFunctional_Elliptic(pb_type(4),pb_type(3),mymesh,...
%                                                           uh,qh,uhat,source_f,...
%                                                           vh,ph,vhat,source_g,...
%                                                           GQ1DRef_pts,GQ1DRef_wts,para.order,para.tau ,0);
%                     Jh_list(ii) = Jh;
%                     Jh_AC_list(ii) = Jh_AC;
%                     ACh_list(ii) = ACh;
%                     % Calculate error
%                     if err_cal_flag==1
%                         % error of uh, qh
%                         [err_uh_list(ii),~] = L2Error_scalar(mymesh,uh,...
%                             GQ1DRef_pts,GQ1DRef_wts,0,...
%                             para.order,uexact);
% 
%                         [err_qh_list(ii),~] = L2Error_vector(mymesh,qh,...
%                             GQ1DRef_pts,GQ1DRef_wts,0,...
%                             para.order,qexact_1,qexact_2);
% 
                        [J_exact,J_exact_elewise_tri]= Exact_Functional(pb_type(4),para.pb_parameters,...
                                    mymesh,GQ1DRef_pts,GQ1DRef_wts);
%                         err_Jh_fine_list(ii) = J_exact - Jh;
%                         err_Jh_AC_fine_list(ii) = J_exact - Jh_AC;
% 
%                         
%                         % error vh, ph;
%                         
%                         [err_vh_list(ii),~] = L2Error_scalar(mymesh,vh,...
%                             GQ1DRef_pts,GQ1DRef_wts,0,...
%                             para.order,vexact);
%                         %PlotElementWiseValue(mymesh,err_vh_elewise,'err-vh elementwise' );
% 
                     end


                end
                % post-processing
                if para.post_process_flag == 1% && Niter>1

                    if strcmp(pb_type(2),'0') && (strcmp(pb_type(4),'1') ||strcmp(pb_type(4),'2'))
                        %%  poission source problem
% 
%                         [uhstar,qhstar] = HDG_Local_Postprocess_Elliptic(mymesh,para.order,uh,qh,uhat,para.tau,GQ1DRef_pts,GQ1DRef_wts);
%                         [vhstar,phstar] = HDG_Local_Postprocess_Elliptic(mymesh,para.order,vh,ph,vhat,para.tau,GQ1DRef_pts,GQ1DRef_wts);
%                         
                        % compute uhstar as HDG k+1
                        
                        [uhstar,~,~] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                        para.order+1, para.tau,source_f,uD,uN);

                        [vhstar,~,~] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                        para.order+1, para.tau,source_g,vD,vN);


                    end
                end
            end

        %----------------------------------------------------------------------
        % step 4: Get uh_proj, vh_proj by L2 projection of uh, vh to TH
        %         Get uH_square, vH_square by L2 projection of uH_triangle,
        %         vH_triangle to TH
        %----------------------------------------------------------------------
            err_Jh_coarse_list(jj) = abs(err_Jh_fine_list(Niter));
            err_Jh_AC_coarse_list(jj) = abs(err_Jh_AC_fine_list(Niter));
            ACh_coarse_list(jj) = ACh_list(Niter);
            
            % ------------------------------------------------------------
            % get basic mesh info
            [GQ_x, GQ_y, hx, hy,Nx_coarse,Ny_coarse] = GetPhyGQPts(para.structure_flag,para.dom_type,...
                    h_coarse,GQ1DRef_pts, para.geo_parameters{:});
                
            if Nx_coarse ~= Ny_coarse
                error('Nx!= Ny but the code assume square domain')
            end
            
            order1_dist = (2*poly_k+1)*hx;
            %order1_dist = 0.3;

            N_corner_x = ceil(order1_dist/hx);
            N_corner_y = ceil(order1_dist/hy);
            
            % ------------------------------------------------------------
            % build an outer mesh
            order1_elements = get_order1_elements(1,1,Nx_coarse,Ny_coarse,N_corner_x,N_corner_y);
            mask = build_mesh_mask(para.dom_type, Nx_coarse, Ny_coarse, N_bd, order1_elements);
            out_flag = 1;
            outer_mesh = build_mesh_by_mask(para.dom_type,coarse_mesh, mask, out_flag, Nx_coarse,Ny_coarse);
            %outer_mesh.PlotElement(0)
            
            outer_area = (outer_mesh.num_elements)/2 *hx*hy;
            % inner area = 1 - outter area
            inner_area = 1-outer_area;
            % ------------------------------------------------------------
            
            
            % ------------------------------------------------------------
            % project uH_triagnle(Pk) to uH_square (Qk)
            %{
            uH_square = GetUhL2ProjectionCoarseMesh(poly_k,coarse_mesh,coarse_mesh,...
                    uH_triangle,GQ1DRef_pts, GQ1DRef_wts,hx,hy);
            vH_square = GetUhL2ProjectionCoarseMesh(poly_k,coarse_mesh,coarse_mesh,...
                    vH_triangle,GQ1DRef_pts, GQ1DRef_wts,hx,hy);
                
        
            % project uh_final to TH
            %uh_proj = GetUhL2ProjectionCoarseMesh(poly_k,coarse_mesh,finer_mesh,...
            %        uh_final,GQ1DRef_pts, GQ1DRef_wts,hx,hy);
            %uh_proj_GQpts = GetUhProjGQpts(uh_proj,poly_proj,GQ1DRef_pts);

            uH_square_GQpts = GetUhProjGQpts(uH_square,poly_proj,GQ1DRef_pts); 
            vH_square_GQpts = GetUhProjGQpts(vH_square,poly_proj,GQ1DRef_pts); 
            %}
            % ------------------------------------------------------------
            if para.post_process_flag == 1
                % Project P_{k+1} tp Q_{k+1}
                uH_local_post_square = GetUhL2ProjectionCoarseMesh(poly_k+1,coarse_mesh,coarse_mesh,...
                        uhstar,GQ1DRef_pts, GQ1DRef_wts,hx,hy);
                vH_local_post_square = GetUhL2ProjectionCoarseMesh(poly_k+1,coarse_mesh,coarse_mesh,...
                        vhstar,GQ1DRef_pts, GQ1DRef_wts,hx,hy); 

                %r_ind = QkHierarchical_indx(poly_k+1, poly_k);
                % extrac Q_{k} as uH_square
                % Nk = k+1
                uH_square = uH_local_post_square(1:poly_k+1,1:poly_k+1,:);
                vH_square  = vH_local_post_square(1:poly_k+1,1:poly_k+1,:);
                uH_square_GQpts = GetUhProjGQpts(uH_square,poly_proj,GQ1DRef_pts); 
                vH_square_GQpts = GetUhProjGQpts(vH_square,poly_proj,GQ1DRef_pts); 

            end
           
           
            if flag_conv == 1
                
                % Convolution Postprocessing
                M_uh_conv = ConvolutionFiltering(para.dom_type,spline_degree,...
                        uH_square, Nx_coarse, Ny_coarse, mask, Conv_matrix);
                M_vh_conv = ConvolutionFiltering(para.dom_type,spline_degree,...
                        vH_square, Nx_coarse, Ny_coarse, mask, Conv_matrix);
                    
                if para.post_process_flag == 1
                    M_uh_local_post_conv = ConvolutionFiltering(para.dom_type,spline_degree,...
                            uH_local_post_square, Nx_coarse, Ny_coarse, mask, Conv_matrix_local_post);
                    M_vh_local_post_conv = ConvolutionFiltering(para.dom_type,spline_degree,...
                            vH_local_post_square, Nx_coarse, Ny_coarse, mask, Conv_matrix_local_post);
                    uH_local_post_star_inner = DG_class_square(2*poly_proj+2,hx,hy,Nx_coarse,Ny_coarse,GQ1DRef_pts,mask);
                    uH_local_post_star_inner = uH_local_post_star_inner.set_dg_gq_pts(M_uh_local_post_conv);

                    vH_local_post_star_inner = DG_class_square(2*poly_proj+2,hx,hy,Nx_coarse,Ny_coarse,GQ1DRef_pts,mask);
                    vH_local_post_star_inner = vH_local_post_star_inner.set_dg_gq_pts(M_vh_local_post_conv);

                    [grad_uH_local_post_star_x_GQ, grad_uH_local_post_star_y_GQ]...
                        = uH_local_post_star_inner.get_DG_grad_pts(GQ1DRef_pts,GQ1DRef_pts);
                    [grad_vH_local_post_star_x_GQ, grad_vH_local_post_star_y_GQ]...
                        = vH_local_post_star_inner.get_DG_grad_pts(GQ1DRef_pts,GQ1DRef_pts);
                end
                    

                %Mh_proj = ConvolutionFiltering(para.dom_type,spline_degree,...
                 %       uh_proj, Nx_coarse, Ny_coarse, mask, Conv_matrix);
                    
                % creat DG solution class
                % 2k+1
                uH_star_inner = DG_class_square(2*poly_proj+1,hx,hy,Nx_coarse,Ny_coarse,GQ1DRef_pts,mask);
                uH_star_inner = uH_star_inner.set_dg_gq_pts(M_uh_conv);
                func_uH_star_inner = @(p) uH_star_inner.eval(p(:,1),p(:,2));
                
                vH_star_inner = DG_class_square(2*poly_proj+1,hx,hy,Nx_coarse,Ny_coarse,GQ1DRef_pts,mask);
                vH_star_inner = vH_star_inner.set_dg_gq_pts(M_vh_conv);
                func_vH_star_inner = @(p) vH_star_inner.eval(p(:,1),p(:,2));
                
                
                
                % Compute Functional based on uhstar and vhstar
                % Square elements in inner and triangle elements in outer
                
                source_f_GQ_pts = GetUexactGQpts(source_f, GQ_x, GQ_y);
                source_g_GQ_pts = GetUexactGQpts(source_g, GQ_x, GQ_y);
                uH_star_GQ = uH_star_inner.dg_gq_pts;
                vH_star_GQ = vH_star_inner.dg_gq_pts;
                
                % Get grad uh in x and y component
                [grad_uH_star_x_GQ, grad_uH_star_y_GQ] ...
                    = uH_star_inner.get_DG_grad_pts(GQ1DRef_pts,GQ1DRef_pts);
                
                % Get grad vh in x and y component
                [grad_vH_star_x_GQ, grad_vH_star_y_GQ]...
                    = vH_star_inner.get_DG_grad_pts(GQ1DRef_pts,GQ1DRef_pts);
                
%                 [Jh_star_inner, ACh_star_inner] = ...
%                     LinearFunctional_uhstar_inner(uH_star_inner, source_f_GQ_pts,...
%                     vH_star_inner,source_g_GQ_pts,...
%                     GQ1DRef_pts,GQ1DRef_wts, hx, hy);
                
                [Jh_star_inner, ACh_star_inner, Jh_inner_elewise, ACh_inner_elewise] = ...
                    LinearFunctional_inner_gradient(...
                    uH_star_GQ,grad_uH_star_x_GQ,grad_uH_star_y_GQ, source_f_GQ_pts,...
                    vH_star_GQ,grad_vH_star_x_GQ,grad_vH_star_y_GQ, source_g_GQ_pts,...
                    uH_star_inner.Nx,uH_star_inner.Ny,uH_star_inner.mask,...
                    GQ1DRef_wts, uH_star_inner.hx, uH_star_inner.hy);
                
                if para.post_process_flag == 1
                    % compute a heuristic error estimator
                    err_inner_estimate = FuncErrInnerEstimate(...
                        uH_star_inner.mask,uH_star_inner.Nx,uH_star_inner.Ny,...
                        uH_star_inner.hx, uH_star_inner.hy,...
                        grad_uH_local_post_star_x_GQ,grad_uH_local_post_star_y_GQ,...
                        grad_uH_star_x_GQ,grad_uH_star_y_GQ,...
                        grad_vH_local_post_star_x_GQ,grad_vH_local_post_star_y_GQ,...
                        grad_vH_star_x_GQ, grad_vH_star_y_GQ,...
                        GQ1DRef_wts...
                        );
                    ave_err_inner_estimate = err_inner_estimate/(inner_area);
                else
                    err_inner_estimate = 0;
                end
                
%                 J_exact_inner = VisulaizeElewiseError(...
%                     J_exact_elewise_tri,Jh_inner_elewise,ACh_inner_elewise,...
%                     uH_star_inner.mask,uH_star_inner.Nx,uH_star_inner.Ny);
%                 
                
                %%%%%%%%%%%%%%%%%%%%%%
                % compute error in the interior domain
%                 qexact_1_GQ = GetUexactGQpts(qexact_1, GQ_x, GQ_y);
%                 qexact_2_GQ = GetUexactGQpts(qexact_2, GQ_x, GQ_y);
%                 pexact_1_GQ = GetUexactGQpts(pexact_1, GQ_x, GQ_y);
%                 pexact_2_GQ = GetUexactGQpts(pexact_2, GQ_x, GQ_y);
%                 [Eh_inner,Eh_inner_elmtwise] = LinearFunctonal_inner_error(...
%                     qexact_1_GQ,qexact_2_GQ,-grad_uH_star_x_GQ,-grad_uH_star_y_GQ,...
%                     pexact_1_GQ,pexact_2_GQ, -grad_vH_star_x_GQ, -grad_vH_star_y_GQ,...
%                     uH_star_inner.Nx,uH_star_inner.Ny,uH_star_inner.mask,...
%                     GQ1DRef_wts, uH_star_inner.hx, uH_star_inner.hy );
%                 fprintf('Eh_inner: %.2e   Eh_inner_esti: %.2e\n',Eh_inner,err_inner_estimate);
%                 
%                 
                %%%%%%%%%%%%%%%%%%%%%%
                % solve adaptively on the outer mesh
                [uh2k_outer,~,~,vh2k_outer,~,~,outer_adap_mesh,Jh_outer,ACh_outer] =...
                    Functional_Outer_Driver(outer_mesh, para,poly_outer_k,...
                    func_uH_star_inner,func_vH_star_inner,...%uexact,func_vH_star_inner,...%
                    N_outer_adap_steps,J_exact- Jh_star_inner,ACh_star_inner,err_inner_estimate,ave_err_inner_estimate,outer_area);
                
                % compute error
                err_Jh_star_list(jj) = abs(J_exact- Jh_star_inner-Jh_outer);
                ACh_star_list(jj) = ACh_star_inner+ACh_outer;
                
                err_Jh_ACh_star_list(jj) = abs(J_exact...
                    - Jh_star_inner - ACh_star_inner...
                    - Jh_outer - ACh_outer);
                
                
                
 
                % some plots
                flag_2D_plot = 0;
                if flag_2D_plot == 1
                    Plot2D(para.dom_type, GQ_x, GQ_y, M_uh_conv, "$u_H^*$ on coarse mesh ",0,"")
                end
                
            end

            if err_cal_flag
                
                % ------------------------------------------------------------
                % set up
                uexact_GQ_pts = GetUexactGQpts(uexact, GQ_x, GQ_y);
                
                % L2 error on the whole domain (square mesh)
                err_uH = L2Error_scalar_Square(uH_square_GQpts,...
                        uexact_GQ_pts, GQ1DRef_wts, hx, hy);
                
                err_uH_list(jj) =   err_uH;
                
                % ------------------------------------------------------------
 
                flag_plot_diff_no_postprocess = 0;
                if flag_plot_diff_no_postprocess == 1 && jj== Ncoarse

                    name_text = "k"+num2str(poly_k)+"_Adapt_Mesh"+num2str(Niter) + "_err_uH";
                    Plot2D(para.dom_type, GQ_x, GQ_y, uH_square_GQpts,... % uexact_GQ_pts-
                        "$u - u_{H,square}$ on coarse mesh ",save_flag,name_text)

                    name_text = "k"+num2str(poly_k)+"_Adapt_Mesh"+num2str(Niter) + "_err_uhproj";
                    Plot2D(para.dom_type, GQ_x, GQ_y, uh_proj_GQpts,... % uexact_GQ_pts-
                        "$u - u_{h,proj}$ on coarse mesh ",save_flag,name_text)
                end 

                
                if flag_conv == 1
                    
                    % L2 error of uhstar on inner domain
                    % now only consider the inner domain
                    uexact_GQ_pts = get_inner_domain_data(para.dom_type, uexact_GQ_pts,Nx_coarse, Ny_coarse, mask);
                    err_uHstar = L2Error_scalar_Square(M_uh_conv,...
                            uexact_GQ_pts, GQ1DRef_wts, hx, hy);
                    err_uH_star_list(jj) = err_uHstar/sqrt(inner_area);
                    
                
                    
                    ndof_inner = GetInnerDof(coarse_mesh, outer_mesh, poly_k);
                    ndof_outer = GetDof(outer_adap_mesh, poly_outer_k);
                    

                    % evaluate the error
                    [err_uh2k_outer,~] = L2Error_scalar(outer_adap_mesh,uh2k_outer,...
                                GQ1DRef_pts,GQ1DRef_wts,0,...
                                poly_outer_k,uexact);
                            
                    % L2 error of uhstar on the whole domain
                    err_uH_star_comp_list(jj) = sqrt(err_uHstar^2 + err_uh2k_outer^2);
                    
                    mesh_comp_list(jj) = ndof_inner + ndof_outer;
                    
                    fprintf('err_inner: %.2e err_outer: %.2e\n',err_uHstar ,err_uh2k_outer);
                    fprintf('err_inner/sqrt(inn_area): %.2e, err_outer/sqrt(out_area): %.2e\n',...
                        err_uHstar/sqrt(inner_area),err_uh2k_outer/sqrt(1-inner_area));
                    fprintf('dof_inner %d,  dof_outer: %d\n', ndof_inner, ndof_outer)

                    flag_plot_diff = 0;
                    if flag_plot_diff == 1 && jj== Ncoarse

                        name_text = "k"+num2str(poly_k)+"_Adapt_Mesh"+num2str(Niter) + "_err_uHstar";

                        title_text = "$u- u_H^*$, L2 error: " + sprintf('%.2e',err_uHstar);

                        Plot2D(para.dom_type, GQ_x, GQ_y, M_uh_conv,... % uexact_GQ_pts-
                            title_text,save_flag,...
                            name_text)

                        name_text = "k"+num2str(poly_k)+"_Adapt_Mesh"+num2str(Niter) + "_err_uhprojstar";

                        title_text = "$u - u_{h,proj}^*$, L2 error: "+ sprintf('%.2e',err_uhstar);
                        Plot2D(para.dom_type, GQ_x, GQ_y, Mh_proj,... % uexact_GQ_pts-
                            title_text,save_flag,...
                            name_text) 

                    end
                end

            end

         
        end
        
        % ------------------------------------------------------------
        % Report result
        % ------------------------------------------------------------
        order_uH        = GetOrderH(h0_list,err_uH_list);
        ReportTable('h', h0_list,...
                    'err_uH_Om', err_uH_list, 'order',order_uH)
                
        if flag_conv == 1
            order_uH_star   = GetOrderH(h0_list,err_uH_star_list);
            ReportTable('h', h0_list,...
                        'err_uH*_inner', err_uH_star_list, 'order',order_uH_star)
 
            order_uH_star_comp = GetOrder(mesh_comp_list,err_uH_star_comp_list);
            
            ReportTable('DoF',mesh_comp_list,...
                'err_uH*_comp_Om',err_uH_star_comp_list,'order',order_uH_star_comp)
            
            
            
            
        end
        
        order_Jh = GetOrderH(h0_list,err_Jh_coarse_list);
        order_Jh_AC = GetOrderH(h0_list,err_Jh_AC_coarse_list);
        ReportTable('h', h0_list,...
                'J-J(uh)',err_Jh_coarse_list,'order',order_Jh,...
                'ACh',ACh_coarse_list, 'ratio', err_Jh_coarse_list./ACh_coarse_list,...
                'J-J(uh)-ACh',err_Jh_AC_coarse_list,'order',order_Jh_AC);
            
        if flag_conv == 1
            order_Jhstar = GetOrder(mesh_comp_list,err_Jh_star_list);
            order_Jhstar_ACh = GetOrder(mesh_comp_list,err_Jh_ACh_star_list);
%             order_Jhstar = GetOrderH(h0_list,err_Jh_star_list);
%             order_Jhstar_ACh = GetOrderH(h0_list,err_Jh_ACh_star_list);
            order_ACh_star = GetOrderH(h0_list,abs(ACh_star_list));
            ReportTable('DoF',mesh_comp_list,...
                'J-J(uh*)',err_Jh_star_list,'order',order_Jhstar,...
                'ACh*',abs(ACh_star_list),'err/Ach*', err_Jh_star_list./abs(ACh_star_list),...%'order', order_ACh_star,...%
                'J-J(uh*)-ACh*',err_Jh_ACh_star_list,'order',order_Jhstar_ACh)
            
            figure;
            plot(0.5*log10(mesh_comp_list),log10(err_Jh_star_list),'--bo');
            hold on
            plot(0.5*log10(mesh_comp_list),log10(abs(ACh_star_list)),'--ks');
            plot(0.5*log10(mesh_comp_list),log10(err_Jh_ACh_star_list),'--rx');
            legend('J-J(uh*)', 'Ach*','J-J(uh*)-ACh*');
        end
            
            
        
    end
 
end