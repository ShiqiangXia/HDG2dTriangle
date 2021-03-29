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
                 [mu,epsilon,omg]...
                     =MyParaParse(para.pb_parameters,'mu','epsilon','omg');
                 
                 [lamh,wh_Neig,uh_Neig,ph_Neig,hat_var_Neig]...
                    = HDG_EigPbSolver_Maxwell(pb_type(3),mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                    para.order,para.tau,para.tau,mu,epsilon,omg,...
                    Neig,Max_iter,Tol_eig);
                
                lamh_list(ii,:) = lamh;
                
                err_lamh_list(ii,:) = EigenError(pb_type(3),lamh,para.dom_type,para.geo_parameters);
                
                if err_cal_flag 
                    % which eigenfunction we want to compute error for. 
                    [tag_eig] = MyParaParse(para.pb_parameters,'tag_eig');
                    wh = wh_Neig(:,:,:,tag_eig);
                    uh = uh_Neig(:,:,:,tag_eig);
                    ph = ph_Neig(:,:,:,tag_eig);
                    hat_var = hat_var_Neig(:,:,:,tag_eig);
                end
                
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
            
            if strcmp(pb_type(2),'1')
                tag_text = 'eig_';
                temp_eig=ParseEigenError(mesh_list,err_lamh_list,tag_text);
                ReportTable('DOF', mesh_list,...
                    temp_eig{:})
            end
            

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
        
     elseif strcmp(pb_type(1),'2')
         
         if strcmp(pb_type(2),'0') % source problem
             
         elseif strcmp(pb_type(2),'1') % eigenproblem
            err_lamh_AC_list = zeros(Niter,Neig,numeric_t);
            lamh_AC_list = zeros(Niter,Neig,numeric_t);
            lamh2_list = zeros(Niter,Neig,numeric_t);
            err_lamh2_list = zeros(Niter,Neig,numeric_t);
            ACh_list = zeros(Niter,Neig,numeric_t);
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
            mymesh.Plot(0);
            
            % Solve -------------------------------------------------------
            if strcmp(pb_type(2),'0') % source problem
            elseif strcmp(pb_type(2),'1')
                % Solve eigen problem
                [mu,epsilon,omg]...
                =MyParaParse(para.pb_parameters,'mu','epsilon','omg');

                [lamh,wh_Neig,uh_Neig,ph_Neig,hat_var_Neig]...
                = HDG_EigPbSolver_Maxwell(pb_type(3),mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                para.order,para.tau,para.tau,mu,epsilon,omg,...
                Neig,Max_iter,Tol_eig);

                lamh_list(ii,:) = lamh;

                err_lamh_list(ii,:) = EigenError(pb_type(3),lamh,para.dom_type,para.geo_parameters);
                % use non-linear functional formula to compute lambdah

                [lamh2,lamh_AC,ACh,ACh_Neig_elewise_list]=...
                EigenvalueFunctional_Maxwell(pb_type(3),mymesh,Neig,...
                wh_Neig,uh_Neig,hat_var_Neig,...
                GQ1DRef_pts,GQ1DRef_wts,para.order,para.tau,para.tau,mu,epsilon,omg ); 

                lamh2_list(ii,:) = lamh2;
                lamh_AC_list(ii,:) = lamh_AC;
                ACh_list(ii,:) = ACh;

                err_lamh_AC_list(ii,:) = EigenError(pb_type(3),lamh_AC,para.dom_type,para.geo_parameters);
                err_lamh2_list(ii,:)  = EigenError(pb_type(3),lamh2,para.dom_type,para.geo_parameters);
                % adaptivity based on which eigenvalue
                [tag_eig] = MyParaParse(para.pb_parameters,'tag_eig');
                
                ACh_elewise_list = ACh_Neig_elewise_list(:,tag_eig);
                
            end
            
            % Posterior error estimate if needed--------------------------- 
            if para.refine_flag > 0
                mark_flag = 0; % 1: bulk marking strategy Dorfler , 0: max marking strategy
                [tol_adp,percent] = MyParaParse(para.extra_parameters,'tol_adp','percent');
                marked_elements = ACh_ErrEstimate(ACh_elewise_list,tol_adp,percent,mark_flag);
                
                % Plot estimator
                if ii <= Niter
                    title_text = append('ACh element-wise, mesh: ',num2str(ii));
                    PlotElementWiseValue(mymesh,abs(ACh_elewise_list),title_text);
                end
                
            end
            
            
        end
        
        if para.report_flag==1
            if strcmp(pb_type(2),'1')

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
                temp_eig{:});
            
                figure;
                plot(0.5*log10(mesh_list),log10(err_lamh2_list(:,tag_eig)),'--bo',...
                    0.5*log10(mesh_list),log10(err_lamh_AC_list(:,tag_eig)),'--kx',...
                    0.5*log10(mesh_list),log10(abs(ACh_list(:,tag_eig))),'--rs');
                legend('Err-lamh','Err-lamh-AC','ACh')
                title('Log plot of errors and estimator for the 1st eigenvalue');
            

            end

        end
            
              
            
             
            
        

     else
         error('Pb type is not incorrect, please double check and see Parameter obj')
     end
    
    
end