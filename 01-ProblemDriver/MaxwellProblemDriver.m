function MaxwellProblemDriver(para)
    % This is the main problem driver for maxwell equations.
    % unknowns: wh,uh,ph, hat_var ( uhat_t, phat)
    
    pb_type = num2str(para.pb_type);
    Niter = para.Niter;
    
    %% step 1. Set varibales to store results %%%%%%%%%%%%%%%%%%%%%
    err_cal_flag = para.err_cal_flag;
    mesh_list = zeros(Niter,1,numeric_t);
    
    if err_cal_flag == 1
        err_wh_list = zeros(Niter,1,numeric_t);
        err_uh_list = zeros(Niter,1,numeric_t);
        err_ph_list = zeros(Niter,1,numeric_t);
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
    
 
     
     %% step 2: Determine what problem we are solving.
     
     if strcmp(pb_type(1),'1')
         % Solve PDE problem
         
         %%%%%%% step 3. Iterative  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cprintf('blue','--------------------------------\n')
        cprintf('blue','Start solving PDE problem\n')
        for ii = 1:Niter
            cprintf('blue','Mesh %d ... \n',ii)
            % build mesh --------------------------------------------------
            if ii == 1 
                % build first mesh
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
                 [source_j_1,source_j_2,uxn_D,mu,epsilon,omg]...
                     =MyParaParse(para.pb_parameters,'source_j_1','source_j_2',...
                     'uxn_D','mu','epsilon','omg');
                 
                 [wh,uh,ph,hat_var] ...
                 = HDG_SourcePbSolver_Maxwell(pb_type(3),mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                        para.order,para.tau,para.tau,...
                        source_j_1,source_j_2,uxn_D,mu,epsilon,omg);
             elseif strcmp(pb_type(2),'1')
                 % Solve eigen problem
             else
                error('pb type not implemented yet')
             end
             
             % Calculate function Error ---------------------------------------------
            
            mesh_list(ii) = GetDof(mymesh, para.order);
            
            if err_cal_flag
                [wexact,uexact_1,uexact_2,pexact]=MyParaParse(para.pb_parameters,'wexact','uexact_1','uexact_2','pexact');
                
                [err_wh_list(ii),err_wh_elewise] = L2Error_scalar(mymesh,wh,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    para.order,wexact);
                
                PlotElementWiseValue(mymesh,err_wh_elewise,'err-wh elementwise' );
                
                [err_uh_list(ii),~] = L2Error_vector(mymesh,uh,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    para.order,uexact_1,uexact_2);
                
                [err_ph_list(ii),~] = L2Error_scalar(mymesh,ph,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    para.order,pexact);
                
            end
            
            % visualization -----------------------------------------------
            if ii == Niter 
                mymesh.Plot(1);
            end
            
        end
            
        % Step 4. Report reulsts%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if para.report_flag==1
            ReportProblem(para)

            if  err_cal_flag == 1
                order_wh = GetOrder(mesh_list,err_wh_list);
                order_uh = GetOrder(mesh_list,err_uh_list);
                order_ph = GetOrder(mesh_list,err_ph_list);
                ReportTable('DOF', mesh_list,...
                    'err_wh',err_wh_list,'order', order_wh,...
                    'err_uh',err_uh_list,'order',order_uh,...
                    'err_ph',err_ph_list,'order', order_ph);

            end

        end
              
            
             
            
        

     else
         error('Pb type is not incorrect, please double check and see Parameter obj')
     end
    
    
end