function EllipticProblemDriver(para)
    % This is the main problem driver for elliptic equations.
    % three unknowns: qh, uh, uhat
    
    %% Step 1 : Set up some varibales
    pb_type = num2str(para.pb_type);
    Niter = para.Niter;
    posterior_estimate_method = 2;
    
    %%%%%% step 1. Set varibales to store results %%%%%%%%%%%%%%%%%%%%%
    err_cal_flag = para.err_cal_flag;
    mesh_list = zeros(Niter,1,numeric_t);
    
    if err_cal_flag == 1
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
    
    %%  step 2: Determine what problem we are solving.
    
    if strcmp(pb_type(1),'1')
        
        %%  Solve PDE problem

        % step 3. Iterative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cprintf('blue','--------------------------------\n')
        cprintf('blue','Start solving PDE problem\n')
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
                       
            % Solve -------------------------------------------------------

            if strcmp(pb_type(2),'0')
                % Solve source problem
                
                [source_f,uD,uN]=MyParaParse(para.pb_parameters,'source_f','uD','uN');
                [uh,qh,uhat] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau,source_f,uD,uN);
                
            % -------------------------------------------------------------

            elseif strcmp(pb_type(2),'1')
                % Solve eigen problem
                [lamh,uh_Neig,qh_Neig,uhat_Neig] = HDG_EigPbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                    para.order,para.tau, Neig,Max_iter,Tol_eig);
                lamh_list(ii,:) = lamh;

                err_lamh_list(ii,:) = EigenError(pb_type(3),lamh,para.dom_type,para.geo_parameters);
          
                if err_cal_flag 
                    % which eigenfunction we want to compute error for. 
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
            
            % Calculate function Error ---------------------------------------------
            
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
            end
            
            % post processed error 
            if para.post_process_flag == 1
                
                [uhstar,qhstar] = HDG_Local_Postprocess_Elliptic(mymesh,para.order,uh,qh,uhat,para.tau,GQ1DRef_pts,GQ1DRef_wts);
                
                if err_cal_flag == 1
                    err_uhstar_list(ii) = L2Error_scalar(mymesh,uhstar,...
                        GQ1DRef_pts,GQ1DRef_wts,1,...
                    para.order,uexact);

                    err_qhstar_list(ii)= L2Error_vector(mymesh,qhstar,...
                       GQ1DRef_pts,GQ1DRef_wts,1,...
                    para.order,qexact_1,qexact_2);
                end
            end
            % ------------------------------------------------------------- 

            % visualization -----------------------------------------------
            if ii == Niter && para.visualize_flag==1
                mymesh.Plot(1);
            end
            if ii == Niter && para.visualize_flag==1
                basis_flag = 0;
                My2DTriPlot(mymesh,uh,para.order, GQ1DRef_pts,basis_flag );
            end
            % ------------------------------------------------------------- 
            
            
        end
        
        % Step 4. Report reulsts%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if para.report_flag==1 
            
            ReportProblem(para) 
            
            % report eigenvalue results
            if strcmp(pb_type(2),'1')
                tag_text = 'eig_';
                temp_eig=ParseEigenError(mesh_list,err_lamh_list,tag_text);
                ReportTable('DOF', mesh_list,...
                    temp_eig{:})
            end
            
            % report uh, qh errors
            if  err_cal_flag == 1
                order_uh = GetOrder(mesh_list,err_uh_list);
                order_qh = GetOrder(mesh_list,err_qh_list);
                ReportTable('DOF', mesh_list,...
                    'err_uh',err_uh_list,'order', order_uh,...
                    'err_qh',err_qh_list,'order',order_qh)

                if para.post_process_flag == 1
                    order_uhstar = GetOrder(mesh_list,err_uhstar_list);
                    order_qhstar = GetOrder(mesh_list,err_qhstar_list);
                    ReportTable('DOF', mesh_list,...
                        'err_uhstar',err_uhstar_list,'order',order_uhstar,...
                        'err_qhstar',err_qhstar_list,'order',order_qhstar )
                end
            end
                 
        end
        % -----------------------------------------------------------------
    
    elseif strcmp(pb_type(1),'2')
        
        %%  Solver Functional problem
        % define some variables to store results
        if strcmp(pb_type(2),'0') % source problem
            Jh_list = zeros(Niter,1,numeric_t);
            err_Jh_list = zeros(Niter,1,numeric_t);
            Jh_AC_list = zeros(Niter,1,numeric_t);
            err_Jh_AC_list = zeros(Niter,1,numeric_t);
            ACh_list = zeros(Niter,1,numeric_t);
            err_vh_list = zeros(Niter,1,numeric_t);
            
            err_terms_sum_list = zeros(Niter,1,numeric_t);
            err_terms_1_list = zeros(Niter,1,numeric_t);
            err_terms_2_list = zeros(Niter,1,numeric_t);
            err_terms_3_list = zeros(Niter,1,numeric_t);
            err_terms_4_list = zeros(Niter,1,numeric_t);
            err_terms_5_list = zeros(Niter,1,numeric_t);
            err_terms_extra_list = zeros(Niter,1,numeric_t);
            
            eterm_1_list = zeros(Niter,1,numeric_t);
            eterm_2_list = zeros(Niter,1,numeric_t);
            eterm_3_list = zeros(Niter,1,numeric_t);
            
            
            est_terms_sum_list = zeros(Niter,1,numeric_t);
            est_terms_1_list = zeros(Niter,1,numeric_t);
            est_terms_2_list = zeros(Niter,1,numeric_t);
            est_terms_3_list = zeros(Niter,1,numeric_t);
            est_terms_4_list = zeros(Niter,1,numeric_t);
            est_terms_5_list = zeros(Niter,1,numeric_t);
            
            post_terms_1_list = zeros(Niter,1,numeric_t);
            post_terms_2_list = zeros(Niter,1,numeric_t);
            post_terms_3_list = zeros(Niter,1,numeric_t);
            post_terms_4_list = zeros(Niter,1,numeric_t);
            
        elseif strcmp(pb_type(2),'1') % eigenproblem
            err_lamh_AC_list = zeros(Niter,Neig,numeric_t);
            lamh_AC_list = zeros(Niter,Neig,numeric_t);
            lamh2_list = zeros(Niter,Neig,numeric_t);
            err_lamh2_list = zeros(Niter,Neig,numeric_t);
            ACh_list = zeros(Niter,Neig,numeric_t);
        else
            error('Wrong problem type')
        end
        % step 3. Iterative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cprintf('blue','--------------------------------\n')
        cprintf('blue','Start solving functional problem\n')

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
            
            
            mymesh.Plot(0); 
            M(ii) = getframe(gcf);
            
            mesh_list(ii) = GetDof(mymesh, para.order);
            
            % Solve -------------------------------------------------------
            if strcmp(pb_type(2),'0') % source problem
                % need to solve Primal and Adjoint two problems
                if strcmp(pb_type(3),'1') % Solve Poission source problem
             
                    [source_f,uD,uN]=MyParaParse(para.pb_parameters,'source_f','uD','uN');
                    [source_g,vD,vN]=MyParaParse(para.pb_parameters,'source_g','vD','vN');
                    
                    [uh,qh,uhat] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau,source_f,uD,uN);
                
                    [vh,ph,vhat] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau,source_g,vD,vN);
                    
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
                
                %PlotElementWiseValue(mymesh,ACh_elewise_list,'ACh-elewise-list');
                

                % Cal funciton Error ------------------------------------------------
                 
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
                    
                    % compute Eh term by term
                    if strcmp(pb_type(4),'1')
                        [err_terms_sum,err_term1,err_term2,err_term3,...
                            err_term4,err_term5,err_term_extra,...
                            eterm1,eterm2,eterm3] ...
                            = Explicit_Functional_Error_Terms(pb_type(4),pb_type(3),para.pb_parameters,...
                            mymesh,...
                            uh,qh,uhat,...
                            vh,ph,vhat,...
                            GQ1DRef_pts,GQ1DRef_wts,para.order,para.tau,0);
                 
                        err_terms_sum_list(ii) = sum(err_terms_sum);
                        err_terms_1_list(ii) =  sum(err_term1);
                        err_terms_2_list(ii) =  sum(err_term2);
                        err_terms_3_list(ii) =  sum(err_term3);
                        err_terms_4_list(ii) =  sum(err_term4);
                        err_terms_5_list(ii) =  sum(err_term5);
                        err_terms_extra_list(ii) =  sum(err_term_extra);
                        
                        
                        eterm_1_list(ii) = sqrt(sum(eterm1));
                        eterm_2_list(ii) = sqrt(sum(eterm2));
                        eterm_3_list(ii) = sqrt(sum(eterm3));
                        
                    
                    end
                    
                    % error vh, ph;
                    [vexact,pexact_1,pexact_2]=MyParaParse(para.pb_parameters,'vexact','pexact_1','pexact_2');
                    [err_vh_list(ii),err_vh_elewise] = L2Error_scalar(mymesh,vh,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order,vexact);
                    %PlotElementWiseValue(mymesh,err_vh_elewise,'err-vh elementwise' );
                    
                end

            elseif strcmp(pb_type(2),'1') % eigenproblem
                %%  only need to solve one eigenvlaue problem
               
                [lamh,uh_Neig,qh_Neig,uhat_Neig] = HDG_EigPbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                para.order,para.tau, Neig,Max_iter,Tol_eig);
            
                lamh_list(ii,:) = lamh;
                err_lamh_list(ii,:) = EigenError(pb_type(3),lamh,para.dom_type,para.geo_parameters);

                % use non-linear functional formula to compute lambdah
                
                [lamh2,lamh_AC,ACh,ACh_Neig_elewise_list]=...
                    EigenvalueFunctional_Elliptic(pb_type(3),mymesh,Neig,...
                    uh_Neig,qh_Neig,uhat_Neig,...
                    GQ1DRef_pts,GQ1DRef_wts,para.order,para.tau ); 
                
                lamh2_list(ii,:) = lamh2;
                lamh_AC_list(ii,:) = lamh_AC;
                ACh_list(ii,:) = ACh;
                
                err_lamh_AC_list(ii,:) = EigenError(pb_type(3),lamh_AC,para.dom_type,para.geo_parameters);
                err_lamh2_list(ii,:)  = EigenError(pb_type(3),lamh2,para.dom_type,para.geo_parameters);
                
                % adaptivity based on which eigenvalue
                [tag_eig] = MyParaParse(para.pb_parameters,'tag_eig');
                
                ACh_elewise_list = ACh_Neig_elewise_list(:,tag_eig);
                
            end
            
            % post-processing
            if para.post_process_flag == 1 
 
                if strcmp(pb_type(2),'0')&& strcmp(pb_type(4),'1') 
                    % poission source problem
                
                    [uhstar,qhstar] = HDG_Local_Postprocess_Elliptic(mymesh,para.order,uh,qh,uhat,para.tau,GQ1DRef_pts,GQ1DRef_wts);
                    [vhstar,phstar] = HDG_Local_Postprocess_Elliptic(mymesh,para.order,vh,ph,vhat,para.tau,GQ1DRef_pts,GQ1DRef_wts);
                    
                    err_uhstar_list(ii) = L2Error_scalar(mymesh,uhstar,...
                        GQ1DRef_pts,GQ1DRef_wts,1,...
                        para.order,uexact);

                    err_qhstar_list(ii)= L2Error_vector(mymesh,qhstar,...
                       GQ1DRef_pts,GQ1DRef_wts,1,...
                        para.order,qexact_1,qexact_2);
                    if posterior_estimate_method == 1
                        %% Eh estimate method 1
                    [est_terms_sum,est_term1,est_term2,est_term3,est_term4,est_term5,...
                        post_term1,post_term2,post_term3,post_term4] ...
                            = Functional_Eh_Estimate_Terms(pb_type(4),pb_type(3),mymesh,...
                            uhstar,qhstar,source_f,...
                            vhstar,phstar,source_g,...
                            uh,qh,uhat,...
                            vh,ph,vhat,...
                            GQ1DRef_pts,GQ1DRef_wts,para.order);
                        
                    est_terms_sum_list(ii) = sum(est_terms_sum);
                    est_terms_1_list(ii) = sum(est_term1);
                    est_terms_2_list(ii) = sum(est_term2);
                    est_terms_3_list(ii) = sum(est_term3);
                    est_terms_4_list(ii) = sum(est_term4);
                    est_terms_5_list(ii) = sum(est_term5);
                    
                    post_terms_1_list(ii) = sqrt(sum(post_term1));
                    post_terms_2_list(ii) = sqrt(sum(post_term2));
                    post_terms_3_list(ii) = sqrt(sum(post_term3));
                    post_terms_4_list(ii) = sqrt(sum(post_term4));
                    
                    elseif posterior_estimate_method == 2
                        %% Eh estimate method 2
                        [est_terms_sum,est1,est2,est3,est4,est5]=...
                            Functional_Eh_Estimate_Residual_method(pb_type(4),pb_type(3),mymesh,...
                            uhstar,qhstar,source_f,...
                            vhstar,phstar,source_g,...
                            uh,qh,uhat,...
                            vh,ph,vhat,...
                            GQ1DRef_pts,GQ1DRef_wts,para.order,para.tau);
                        
                        est_terms_sum_list(ii) = sum(est_terms_sum);
                        est_terms_1_list(ii) = sum(est1);
                        est_terms_2_list(ii) = sum(est2);
                        est_terms_3_list(ii) = sum(est3);
                        est_terms_4_list(ii) = sum(est4);
                        est_terms_5_list(ii) = sum(est5);
                        
                        
                    end
                end            
            end
            
            
            % Posterior error estimate if needed--------------------------- 
            if para.refine_flag > 0
                mark_flag = 0; % 1: bulk marking strategy Dorfler , 0: max marking strategy
                
                if para.post_process_flag == 1
                    estimator_functinal = est_terms_sum+ACh_elewise_list;%+est_terms_sum ;%+%;
                else
                    estimator_functinal = ACh_elewise_list;
                    
                end
                
                [tol_adp,percent] = MyParaParse(para.extra_parameters,'tol_adp','percent');
                marked_elements = ACh_ErrEstimate(estimator_functinal,tol_adp,percent,mark_flag);
                
                % Plot estimator
                if ii <= Niter
%                     title_text = append('ACh element-wise, mesh: ',num2str(ii));
%                     PlotElementWiseValue(mymesh,ACh_elewise_list,title_text,...
%                         est_terms_sum,'Dh',err_terms_sum,'Error Eh element-wise');
%                     
%                     title_text = append('ACh, mesh: ',num2str(ii));
%                     PlotElementWiseValue(mymesh,ACh_elewise_list,title_text);
%                     title_text = append('Dh, mesh: ',num2str(ii));
%                     PlotElementWiseValue(mymesh,est_terms_sum,title_text);
%                     
%                      title_text = append('ACh+Dh, mesh: ',num2str(ii));
%                     PlotElementWiseValue(mymesh,ACh_elewise_list+est_terms_sum,title_text);
%                     
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
        
        if para.visualize_flag==1
            implay(M,2)
        end
        
        % Report reulsts---------------------------------------------------
        
        if para.report_flag==1 
            
            ReportProblem(para) 
            
            if strcmp(pb_type(2),'0')           
                order_Jh = GetOrder(mesh_list,err_Jh_list);
                order_Jh_AC = GetOrder(mesh_list,err_Jh_AC_list);
                
                err_Jh_AC_Dh = err_Jh_AC_list - est_terms_sum_list;
                order_Jh_AC_Dh = GetOrder(mesh_list,err_Jh_AC_Dh);
                
                if para.post_process_flag == 1
                    
                    estimate_err_Jh = ACh_list+est_terms_sum_list;

                    ReportTable('DOF', mesh_list,...
                    'err_Jh',err_Jh_list,'order',order_Jh, ...
                    'ACh+Dh',estimate_err_Jh,'|est/err|', abs(estimate_err_Jh./err_Jh_list),...
                    'err_Jh_AC_Dh',err_Jh_AC_Dh,'order',order_Jh_AC_Dh );

                    ReportTable('DOF', mesh_list,...
                        'err_Jh',err_Jh_list,'order',order_Jh, ...
                        'ACh',ACh_list,'ACh/err', ACh_list./err_Jh_list,...
                        'err_Jh_AC',err_Jh_AC_list,'order',order_Jh_AC,...
                        'Dh',est_terms_sum_list,'Dh/err_Jh_AC',est_terms_sum_list./err_Jh_AC_list);
                else
                    estimate_err_Jh = ACh_list;
                    
                    ReportTable('DOF', mesh_list,...
                        'err_Jh',err_Jh_list,'order',order_Jh, ...
                        'ACh',ACh_list,'ACh/err', ACh_list./err_Jh_list,...
                        'err_Jh_AC',err_Jh_AC_list,'order',order_Jh_AC);
                    
                end

                
                    
                
                %----------------------------------------------------------
                % error terms 
                if strcmp(pb_type(4),'1')
                    
                    %% report explicit Eh terms
                    order_err_terms_sum = GetOrder(mesh_list,err_terms_sum_list);
                    order_err_term_1 = GetOrder(mesh_list,err_terms_1_list);
                    order_err_term_2 = GetOrder(mesh_list,err_terms_2_list);
                    order_err_term_3 = GetOrder(mesh_list,err_terms_3_list);
                    order_err_term_4 = GetOrder(mesh_list,err_terms_4_list);
                    order_err_term_5 = GetOrder(mesh_list,err_terms_5_list);
                    order_err_term_extra = GetOrder(mesh_list,err_terms_extra_list);

                    fprintf('\n');
                    fprintf('Eh = sum err_i\n')
                    fprintf('err_1 = (q-qh,p-ph)\n');
                    fprintf('err_2 = (q-qh,ph+grad_vh) \n');
                    fprintf('err_3 = (qh+grad_uh,p-ph)\n')
                    fprintf('err_4 = <(qhat-q)*n,vh-vhat> \n');
                    fprintf('err_5 =  <uh-uhat,(phat-p)*n>\n');
                    fprintf('err_extra =  <u-uhat,p*n>\n');

                    ReportTable('DOF', mesh_list,...
                        'err_sum',err_terms_sum_list,'order',order_err_terms_sum,...
                        'err_1',err_terms_1_list,'order',order_err_term_1, ...
                        'err_2',err_terms_2_list,'order',order_err_term_2,...
                        'err_3',err_terms_3_list,'order',order_err_term_3);
                    ReportTable('DOF', mesh_list,...
                        'err_4',err_terms_4_list,'order',order_err_term_4,...
                        'err_5',err_terms_5_list,'order',order_err_term_5,...
                        'err_extra',err_terms_extra_list,'order',order_err_term_extra);
                    
                    
                    order_2_and_4 = GetOrder(mesh_list,err_terms_1_list+err_terms_2_list+err_terms_4_list);
                    order_3_and_5_extra = GetOrder(mesh_list,err_terms_1_list+err_terms_3_list+err_terms_5_list+err_terms_extra_list);
                    
                    ReportTable('DOF', mesh_list,...
                        'err_1+err_2+err_4',err_terms_1_list+err_terms_2_list+err_terms_4_list,'order',order_2_and_4,...
                        'err_1+err_3+err_5+err_extra',err_terms_1_list+err_terms_3_list+err_terms_5_list+err_terms_extra_list,'order',order_3_and_5_extra)


%                     ReportTable('(err_2+err_4)/ACh', (err_terms_2_list+err_terms_4_list)./ACh_list,...
%                         '(err_3+err_5+err_extra)/ACh',(err_terms_3_list+err_terms_5_list+err_terms_extra_list)./ACh_list);

%                     
%                     
%                     order_eterm_1 =  GetOrder(mesh_list,eterm_1_list);
%                     order_eterm_2 =  GetOrder(mesh_list,eterm_2_list);
%                     order_eterm_3 =  GetOrder(mesh_list,eterm_3_list);
%                     
%                     fprintf('\n');
%                     fprintf('eterm1 = ||qh+grad_uh||\n');
%                     fprintf('eterm2 = ||(qhat-q)*n||_F\n');
%                     fprintf('eterm3 = ||uhat -uh||_F\n')
%                     ReportTable('DOF', mesh_list,...
%                         'eterm1',eterm_1_list,'order',order_eterm_1,...
%                         'eterm2',eterm_2_list,'order',order_eterm_2,...
%                         'eterm3',eterm_3_list,'order',order_eterm_3);

                    %% report post-processed results
                    if para.post_process_flag == 1
                        
                        if posterior_estimate_method == 1
                            %% 
                            order_pterm_1 = GetOrder(mesh_list,post_terms_1_list);
                            order_pterm_2 = GetOrder(mesh_list,post_terms_2_list);
                            order_pterm_3 = GetOrder(mesh_list,post_terms_3_list);
                            order_pterm_4 = GetOrder(mesh_list,post_terms_4_list);
                            fprintf('\n');
                            fprintf('pterm1 = ||uh*-uh||\n');
                            fprintf('pterm2 = ||qh*-qh||\n');
                            fprintf('pterm3 = ||-graduh* - qh||\n')
                            fprintf('pterm4 = ||-graduh* - qh*||\n')
                            ReportTable('DOF', mesh_list,...
                                'pterm1',post_terms_1_list,'order',order_pterm_1,...
                                'pterm2',post_terms_2_list,'order',order_pterm_2,...
                                'pterm3',post_terms_3_list,'order',order_pterm_3,...
                                'pterm4',post_terms_4_list,'order',order_pterm_4);

                            order_est_terms_sum = GetOrder(mesh_list,est_terms_sum_list);
                            order_est_term_1 = GetOrder(mesh_list,est_terms_1_list);
                            order_est_term_2 = GetOrder(mesh_list,est_terms_2_list);
                            order_est_term_3 = GetOrder(mesh_list,est_terms_3_list);
                            order_est_term_4 = GetOrder(mesh_list,est_terms_4_list);
                            order_est_term_5 = GetOrder(mesh_list,est_terms_5_list);

                            fprintf('\n');
                            fprintf('Dh = sum est_i\n')
                            fprintf('est_1 = (qh*-qh,ph*-ph)\n');
                            fprintf('est_2 = (qh*-qh,ph+grad_vh) \n');
                            fprintf('est_3 = (qh+grad_uh,ph*-ph)\n')
                            fprintf('est_4 = (f-Proj_f, vh* - vh) \n');
                            fprintf('est_5 =  (uh*-uh, g - Proj_g)\n');

                            ReportTable('DOF', mesh_list,...
                                'est_sum',est_terms_sum_list,'order',order_est_terms_sum,...
                                'est_1',est_terms_1_list,'order',order_est_term_1, ...
                                'est_2',est_terms_2_list,'order',order_est_term_2,...
                                'est_3',est_terms_3_list,'order',order_est_term_3);
                             ReportTable('DOF', mesh_list,...
                                'est_4',est_terms_4_list,'order',order_est_term_4,...
                                'est_5',est_terms_5_list,'order',order_est_term_5);

                             fprintf('\n');
                        elseif posterior_estimate_method == 2
                            %%
                            order_est_terms_sum = GetOrder(mesh_list,est_terms_sum_list);
                            order_est_term_1 = GetOrder(mesh_list,est_terms_1_list);
                            order_est_term_2 = GetOrder(mesh_list,est_terms_2_list);
                            order_est_term_3 = GetOrder(mesh_list,est_terms_3_list);
                            order_est_term_4 = GetOrder(mesh_list,est_terms_4_list);
                            order_est_term_5 = GetOrder(mesh_list,est_terms_5_list);
                            
                            fprintf('\n');
                            fprintf('Dh = sum est_i\n')
                            fprintf('est_1 = (f-div.qh,vh*-vh)\n');
                            fprintf('est_2 = (uh*-uh,g-div.ph)\n');
                            fprintf('est_3 = <(qhat-qh)*n, vh-vh*>\n')
                            fprintf('est_4 = <uh-uh*, (phat-ph)*n>\n')
                            fprintf('est_5 = -(graduh*+qh,gradvh*+ph) \n');

                            order_1_3 = GetOrder(mesh_list,est_terms_1_list+est_terms_3_list);
                            order_2_4 = GetOrder(mesh_list,est_terms_2_list+est_terms_4_list);
                            
                            ReportTable('DOF', mesh_list,...
                                 'est_1+est_3',est_terms_1_list+est_terms_3_list,'order',order_1_3,...
                                'est_2+est_4',est_terms_2_list+est_terms_4_list,'order',order_2_4,...
                                'est_5',est_terms_5_list,'order',order_est_term_5);
                            
                            ReportTable('DOF', mesh_list,...
                                'est_sum',est_terms_sum_list,'order',order_est_terms_sum,...
                                'est_1',est_terms_1_list,'order',order_est_term_1, ...
                                'est_2',est_terms_2_list,'order',order_est_term_2);
                             ReportTable('DOF', mesh_list,...
                                 'est_3',est_terms_3_list,'order',order_est_term_3,...
                                'est_4',est_terms_4_list,'order',order_est_term_4);
                            
                            
                            
                        end
                        
                    end                  
          
                end 
                
                %----------------------------------------------------------
                %% Plot log-error 
                figure;            
                if para.post_process_flag == 0
                    plot(0.5*log10(mesh_list),log10(abs(err_Jh_list)),'--bo',...
                            0.5*log10(mesh_list),log10(abs(estimate_err_Jh)),'--rs');
                    legend('Err-Jh','ACh')
                    title('Log plot of errors and estimator');
                    
                elseif para.post_process_flag == 1

                    plot(0.5*log10(mesh_list),log10(abs(err_Jh_list)),'--bo',...
                            0.5*log10(mesh_list),log10(abs(estimate_err_Jh)),'--rs');
                    legend('Err-Jh','ACh+Dh')
                    title('Log plot of errors and estimator');
                  
                end
                %%
                
                if err_cal_flag==1
                    %% report uh,qh error results
                    
                    order_uh = GetOrder(mesh_list,err_uh_list);
                    order_qh = GetOrder(mesh_list,err_qh_list);
                    order_vh = GetOrder(mesh_list,err_vh_list);
                    fprintf('\n\n');
                    fprintf('Report uh,qh error results \n');
                    ReportTable('DOF', mesh_list,...
                        'err_uh',err_uh_list,'order', order_uh,...
                        'err_qh',err_qh_list,'order',order_qh,...
                        'err_vh',err_vh_list,'order', order_vh)
                    
                    if para.post_process_flag == 1
                        order_uhstar = GetOrder(mesh_list,err_uhstar_list);
                        order_qhstar = GetOrder(mesh_list,err_qhstar_list);
                        ReportTable('DOF', mesh_list,...
                            'err_uhstar',err_uhstar_list,'order',order_uhstar,...
                            'err_qhstar',err_qhstar_list,'order',order_qhstar )
                    end
                end
           
                
            elseif strcmp(pb_type(2),'1')
                %% report eigenvalue results
                tag_text = 'eh_HDG_';
                temp_eig=ParseEigenError(mesh_list,err_lamh_list,tag_text);
                
                ReportTable('dof', mesh_list,...
                    temp_eig{:})
                
                tag_text = 'eh_func_';
                temp_eig=ParseEigenError(mesh_list,err_lamh2_list,tag_text);
                ReportTable('dof', mesh_list,...
                    temp_eig{:})
                
                tag_text = 'ACh_eig_';
                temp_eig=ParseEigenError(mesh_list,ACh_list,tag_text);
                ReportTable('dof', mesh_list,...
                    temp_eig{:})
                
                tag_text = 'eh_ac_';
                temp_eig=ParseEigenError(mesh_list,err_lamh_AC_list,tag_text);
                ReportTable('dof', mesh_list,...
                    temp_eig{:})
                
                tag_text = 'ratio_';
                temp_eig=ParseEigenError(mesh_list,err_lamh2_list./ACh_list,tag_text,0);
                fprintf('ration =  err_lambdah/ACh\n');
                ReportTable('dof', mesh_list,...
                    temp_eig{:})
                
                figure;
                plot(0.5*log10(mesh_list),log10(err_lamh2_list(:,tag_eig)),'--bo',...
                    0.5*log10(mesh_list),log10(err_lamh_AC_list(:,tag_eig)),'--kx',...
                    0.5*log10(mesh_list),log10(abs(ACh_list(:,tag_eig))),'--rs');
                legend('Err-lamh','Err-lamh-AC','ACh')
                title('Log plot of errors and estimator');

                
            else
                
            end
        end
               
        
    else
        error('Pb type is not incorrect, please double check and see Parameter obj')
    end
    
end