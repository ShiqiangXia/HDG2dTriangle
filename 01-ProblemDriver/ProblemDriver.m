function ProblemDriver(para)
    % This is the main problem driver.
    
    % step 1: Determine what problem we are solving.
    pb_type = num2str(para.pb_type);
    if strcmp(pb_type(1),'1')
        
        % Solve PDE problem
        Niter = para.Niter;
        
        %%%%%% step 1. Set varibales to store results %%%%%%%%%%%%%%%%%%%%%
        mesh_list = zeros(Niter,1,numeric_t);
        err_uh_list = zeros(Niter,1,numeric_t);
        err_qh_list = zeros(Niter,1,numeric_t);
        if para.post_process_flag == 1
            err_uhstar_list = zeros(Niter,1,numeric_t);
            err_qhstar_list = zeros(Niter,1,numeric_t);
        end
        
        if strcmp(pb_type(2),'1')
            % eigenvalue problem
            % maybe use matlab inputParser later
            extra_para = para.extra_parameters;
            extra_para = reshape(extra_para,[],2)';
            Neig = extra_para{find(strcmp(extra_para,'Neig')),2};
            err_lam_list = zeros(Niter,Neig,numeric_t);
        end
        
        [GQ1DRef_pts,GQ1DRef_wts] = GaussQuad(para.GQ_deg);
        
        %%%%%%% step 2. Iterative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            % -------------------------------------------------------------
            
            % Solve -------------------------------------------------------
            
            % Solve Poission source problem
            if strcmp(pb_type(3),'1') && strcmp(pb_type(2),'0')
                
                [uh,qh,uhat] = HDG_Poission(mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                    para.order, para.tau, para.pb_parameters);
                
            % -------------------------------------------------------------
            % Solve Poission eigen problem
            elseif  strcmp(pb_type(3),'1') && strcmp(pb_type(2),'1')
                
                [lamh,uh,qh,uhat] = HDG_PoissionEig(mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                    para.order, para.tau,para.extra_parameters);
                err_lam_list(ii) = EigenError(para,lamh);
            % -------------------------------------------------------------
            else
                error('pb type not implemented yet')
            end 
            
            % -------------------------------------------------------------
            % -------------------------------------------------------------
            
            % Error -------------------------------------------------------
            mesh_list(ii) = GetDof(mymesh, para.order);
            err_uh_list(ii) = L2Error_scalar(mymesh,uh,...
                GQ1DRef_pts,GQ1DRef_wts,0,...
                para.order,para.pb_parameters);

            err_qh_list(ii) = L2Error_vector(mymesh,qh,...
                GQ1DRef_pts,GQ1DRef_wts,0,...
                para.order,para.pb_parameters);

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
            % -------------------------------------------------------------
            
            % visualization -----------------------------------------------
            if ii == Niter && para.visualize_flag==1
                basis_flag = 0;
                Plot(mymesh,uh,para.order, GQ1DRef_pts,basis_flag );
            end
            % ------------------------------------------------------------- 
            % -------------------------------------------------------------
            
       
        end
        
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
                    'lamh',err_lam_list)
            end
            
        end
        
    
    elseif strcmp(pb_type(1),'2')
        % Solver Functional problem
        
    else
        error('Pb type is not incorrect, please double check and see Parameter obj')
    end
    
end