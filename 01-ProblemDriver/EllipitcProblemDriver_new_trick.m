function EllipitcProblemDriver_new_trick(para)
    
    % This is the main problem driver for elliptic equations.
    % three unknowns: qh, uh, uhat
    
    pb_type = num2str(para.pb_type);
    Niter = para.Niter;
    
    %%%%%% step 1. Set varibales to store results %%%%%%%%%%%%%%%%%%%%%
    err_cal_flag = para.err_cal_flag;
    mesh_list = zeros(Niter,1,numeric_t);
    
    if err_cal_flag == 1
        err_uh_list = zeros(Niter,1,numeric_t);
        err_qh_list = zeros(Niter,1,numeric_t);
        err_uhstar_list = zeros(Niter,1,numeric_t);
        err_qhstar_list = zeros(Niter,1,numeric_t);
    end
    
    [GQ1DRef_pts,GQ1DRef_wts] = GaussQuad(para.GQ_deg);
    
    % step 2: Determine what problem we are solving.
    
    if strcmp(pb_type(1),'1')
        
    elseif strcmp(pb_type(1),'2')
        
        % Solver Functional problem
        if strcmp(pb_type(2),'0') % source problem
            Jh_list = zeros(Niter,1,numeric_t);
            err_Jh_list = zeros(Niter,1,numeric_t);
            Jh_AC_list = zeros(Niter,1,numeric_t);
            err_Jh_AC_list = zeros(Niter,1,numeric_t);
            ACh_list = zeros(Niter,1,numeric_t);
            err_vh_list = zeros(Niter,1,numeric_t);
            err_ph_list = zeros(Niter,1,numeric_t);
            err_vhstar_list = zeros(Niter,1,numeric_t);
            err_phstar_list = zeros(Niter,1,numeric_t);
            
            err_terms_sum_list = zeros(Niter,1,numeric_t);
            err_terms_1_list = zeros(Niter,1,numeric_t);
            err_terms_2_list = zeros(Niter,1,numeric_t);
            err_terms_3_list = zeros(Niter,1,numeric_t);
            err_terms_4_list = zeros(Niter,1,numeric_t);
            err_terms_5_list = zeros(Niter,1,numeric_t);
            err_terms_extra_list = zeros(Niter,1,numeric_t);
            
            est_terms_sum_list = zeros(Niter,1,numeric_t);
            est_terms_1_list = zeros(Niter,1,numeric_t);
            est_terms_2_list = zeros(Niter,1,numeric_t);
            est_terms_3_list = zeros(Niter,1,numeric_t);
            est_terms_4_list = zeros(Niter,1,numeric_t);
            est_terms_5_list = zeros(Niter,1,numeric_t);
            
        elseif strcmp(pb_type(2),'1') % eigenproblem
           
        else
            error('Wrong problem type')
        end
        
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
            
            mesh_list(ii) = GetDof(mymesh, para.order);
            
            % Solve -------------------------------------------------------
            if strcmp(pb_type(2),'0') % source problem
                % need to solve Primal and Adjoint two problems
                
                if strcmp(pb_type(3),'1') % Solve Poission source problem
                    [source_f,uD,uN]=MyParaParse(para.pb_parameters,'source_f','uD','uN');
                    [source_g,vD,vN]=MyParaParse(para.pb_parameters,'source_g','vD','vN');
                    
                    [uhstar,qhstar,uhatstar] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order + 1, para.tau,source_f,uD,uN);
                
                    [vhstar,phstar,vhatstar] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order + 1, para.tau,source_g,vD,vN);
                
                    
                    [uh,qh,uhat] = HDG_ExtractPart_Solutions(para.order + 1,para.order,uhstar,qhstar,uhatstar);
                    [vh,ph,vhat] = HDG_ExtractPart_Solutions(para.order + 1,para.order,vhstar,phstar,vhatstar);
                    
                
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
                
                % Cal funciton Error ------------------------------------------------
                 
                if err_cal_flag==1
                    
                    %  error of u -----------------------------------------
                    [uexact,qexact_1,qexact_2]=MyParaParse(para.pb_parameters,'uexact','qexact_1','qexact_2');

                    [err_uh_list(ii),~] = L2Error_scalar(mymesh,uh,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order,uexact);

                    [err_qh_list(ii),~] = L2Error_vector(mymesh,qh,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order,qexact_1,qexact_2);
                    
                    [err_uhstar_list(ii),~] = L2Error_scalar(mymesh,uhstar,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order+1,uexact);

                    [err_qhstar_list(ii),~] = L2Error_vector(mymesh,qhstar,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order+1,qexact_1,qexact_2);
                    
                    % error of v ------------------------------------------
                    [vexact,pexact_1,pexact_2]=MyParaParse(para.pb_parameters,'vexact','pexact_1','pexact_2');
                    
                    [err_vh_list(ii),err_vh_elewise] = L2Error_scalar(mymesh,vh,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order,vexact);
                    [err_ph_list(ii),~] = L2Error_vector(mymesh,ph,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order,pexact_1,pexact_2);
                    
                    [err_vhstar_list(ii),err_vh_elewise] = L2Error_scalar(mymesh,vhstar,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order+1,vexact);
                    [err_phstar_list(ii),~] = L2Error_vector(mymesh,phstar,...
                        GQ1DRef_pts,GQ1DRef_wts,0,...
                        para.order+1,pexact_1,pexact_2);
                    
                    
                    
                    % error of the functional -----------------------------
                    [err_Jh_list(ii),err_Jh_AC_list(ii),err_Jh_elewise] ...
                        = Error_Functional(pb_type(4),para.pb_parameters,mymesh,GQ1DRef_pts,GQ1DRef_wts,Jh,Jh_AC,Jh_elewise_list);
                    % explicit error terms
                    if strcmp(pb_type(4),'1')
                        
                        [err_terms_sum,err_term1,err_term2,err_term3,err_term4,err_term5,err_term_extra] ...
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
                    
                    end
 
                end
            end
            
            [est_terms_sum,est_term1,est_term2,est_term3,est_term4,est_term5] ...
                            = Functional_Eh_Estimate_Terms_By_Extraction(pb_type(4),pb_type(3),mymesh,...
                            uhstar,qhstar,uhatstar,...
                            vhstar,phstar,vhatstar,...
                            uh,qh,uhat,...
                            vh,ph,vhat,...
                            GQ1DRef_pts,GQ1DRef_wts,para.order);
                        
            est_terms_sum_list(ii) = sum(est_terms_sum);
            est_terms_1_list(ii) = sum(est_term1);
            est_terms_2_list(ii) = sum(est_term2);
            est_terms_3_list(ii) = sum(est_term3);
            est_terms_4_list(ii) = sum(est_term4);
            est_terms_5_list(ii) = sum(est_term5);
            
            % Posterior error estimate if needed--------------------------- 
            if para.refine_flag > 0
                mark_flag = 0; % 1: bulk marking strategy Dorfler , 0: max marking strategy
                estimator_functinal = est_terms_sum+ACh_elewise_list;%+est_terms_sum ;%+%;
                [tol_adp,percent] = MyParaParse(para.extra_parameters,'tol_adp','percent');
                marked_elements = ACh_ErrEstimate(estimator_functinal,tol_adp,percent,mark_flag);

                % Plot estimator
                if ii <= Niter
%                     title_text = append('ACh+Dh, mesh: ',num2str(ii));
%                     PlotElementWiseValue(mymesh,ACh_elewise_list+est_terms_sum,title_text);
%                     
                end
            end

        end
        
        % Report reulsts---------------------------------------------------
        
        if para.report_flag==1 
            ReportProblem(para)
            if strcmp(pb_type(2),'0')
                order_Jh = GetOrder(mesh_list,err_Jh_list);
                order_Jh_AC = GetOrder(mesh_list,err_Jh_AC_list);

                err_Jh_AC_Dh = err_Jh_AC_list - est_terms_sum_list;

                order_Jh_AC_Dh = GetOrder(mesh_list,err_Jh_AC_Dh);

                estimate_err_Jh = ACh_list+est_terms_sum_list;

                ReportTable('DOF', mesh_list,...
                'err_Jh',err_Jh_list,'order',order_Jh, ...
                'ACh+Dh',estimate_err_Jh,'|est/err|', abs(estimate_err_Jh./err_Jh_list),...
                'err_Jh_AC_Dh',err_Jh_AC_Dh,'order',order_Jh_AC_Dh );

                ReportTable('DOF', mesh_list,...
                    'err_Jh',err_Jh_list,'order',order_Jh, ...
                    'ACh',ACh_list,'ACh/err', ACh_list./err_Jh_list,...
                    'Dh',est_terms_sum_list,'Dh/err',est_terms_sum_list./err_Jh_list);
            end
            
            if err_cal_flag==1
                
                    order_uh = GetOrder(mesh_list,err_uh_list);
                    order_qh = GetOrder(mesh_list,err_qh_list);
                    order_vh = GetOrder(mesh_list,err_vh_list);
                    order_ph = GetOrder(mesh_list,err_ph_list);
                    ReportTable('DOF', mesh_list,...
                        'err_uh',err_uh_list,'order', order_uh,...
                        'err_qh',err_qh_list,'order',order_qh,...
                        'err_vh',err_vh_list,'order', order_vh,...
                        'err_ph',err_ph_list,'order', order_ph)
                    
                  
                    order_uhstar = GetOrder(mesh_list,err_uhstar_list);
                    order_qhstar = GetOrder(mesh_list,err_qhstar_list);
                    order_vhstar = GetOrder(mesh_list,err_vhstar_list);
                    order_phstar = GetOrder(mesh_list,err_phstar_list);
                    
                    ReportTable('DOF', mesh_list,...
                        'err_uhstar',err_uhstar_list,'order',order_uhstar,...
                        'err_qhstar',err_qhstar_list,'order',order_qhstar,...
                        'err_vhstar',err_vhstar_list,'order',order_vhstar,...
                        'err_phstar',err_phstar_list,'order',order_phstar)     

            end
        end

        
    end
    
    
    
    
end