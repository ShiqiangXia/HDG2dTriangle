
refine_flag = 1;
Niter_max = 8;
%TOL_list = [1e-5,1e-6,1e-7];
TOL_list = [1e-14];
N_TOL = size(TOL_list,2);
flag_func1 = 0;
flag_func2 = 1;
flag_func3 = 0;

flag_k1    = 0;
flag_k2    = 0;
flag_k3    = 1;

if refine_flag == 0
    fprintf('Uniform Refinement')
elseif refine_flag == 1
    fprintf('Adaptive Refinement\n')
end



  %% functional 1
  if flag_func1>0
    h0 = 0.25;
    
    fprintf('--------------------------------------------------------------\n')
    fprintf('--------------------------------------------------------------\n')
    fprintf('Initial mesh h0 = %f\n',h0);
    fprintf('\n\nTest Linear functinals J(u) = (u,g)\n\n')
    fprintf('\n\nCase: u corner singularity and g=1 smooth\n')

    pb_type = 2011;
    dom_ype = 'Rec';

    for ii = 1:N_TOL
        fprintf('TOL : %.2e\n',TOL_list(ii))
        if flag_k1>0
            fprintf('--------- k = 1 --------------\n')
            main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
            'pb_type',pb_type,'dom_type',dom_ype,'primal',1,'adjoint',1,...
            'post_process_flag',1,'err_cal_flag',1,'tol_adp',TOL_list(ii));
        end
        
        
        if flag_k2>0
            fprintf('--------- k = 2 --------------\n')
            main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
            'pb_type',pb_type,'dom_type',dom_ype,'primal',1,'adjoint',1,...
            'post_process_flag',1,'err_cal_flag',1,'tol_adp',TOL_list(ii));
        end
        
        if flag_k3>0
            fprintf('--------- k = 3 --------------\n')
            main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
            'pb_type',pb_type,'dom_type',dom_ype,'primal',1,'adjoint',1,...
            'post_process_flag',1,'err_cal_flag',1,'tol_adp',TOL_list(ii));

        end
        
    end

  end

  %% functional 2
if flag_func2 >0
    fprintf('--------------------------------------------------------------\n')
    fprintf('--------------------------------------------------------------\n')
    h0 = 0.2;
    fprintf('Initial mesh h0 = %f\n',h0);
    fprintf('\n\nTest Linear functinals J(u) = <q*n, psi>\n\n')
    fprintf('\n\nCase: u corner singularity and v discontinuous \n')
    
    pb_type = 2012;
    dom_ype = 'Rec';

    for ii = 1:N_TOL
        fprintf('TOL : %.2e\n',TOL_list(ii))
        if flag_k1>0
            fprintf('--------- k = 1 --------------\n')
            main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
            'pb_type',pb_type,'dom_type',dom_ype,'primal',1,'adjoint',2,...
            'post_process_flag',1,'err_cal_flag',1,'tol_adp',TOL_list(ii));
        end
        
        
        if flag_k2>0
            fprintf('--------- k = 2 --------------\n')
            main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
            'pb_type',pb_type,'dom_type',dom_ype,'primal',1,'adjoint',2,...
            'post_process_flag',1,'err_cal_flag',1,'tol_adp',TOL_list(ii));
        end
        
            
        if flag_k3>0
            fprintf('--------- k = 3 --------------\n')
            main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
            'pb_type',pb_type,'dom_type',dom_ype,'primal',1,'adjoint',2,...
            'post_process_flag',1,'err_cal_flag',1,'tol_adp',TOL_list(ii));

        end
        


    end
end

 %% functional 3
 if flag_func3>0
    fprintf('--------------------------------------------------------------\n')
    fprintf('--------------------------------------------------------------\n')

    fprintf('\n\nTest Eigenvalue problem\n\n')
    fprintf('Initial mesh h0 = %f, tag_eig = 1;\n',h0);
    fprintf('\n\nCase:  unit square\n')
    h0 = 0.25;
    pb_type = 2110;
    dom_ype = 'L';
    Niter_max = 7;
    %TOL_list = [5*1e-3,5*1e-4,5*1e-5];
    N_TOL = size(TOL_list,2);

    for ii = 1:N_TOL
        fprintf('TOL : %.2e\n',TOL_list(ii))

        if flag_k1>0
            fprintf('--------- k = 1 --------------\n')
            main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
            'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
            'post_process_flag',1,'err_cal_flag',1,'tol_adp',TOL_list(ii));
        end
        

        if flag_k2>0
            fprintf('--------- k = 2 --------------\n')
            main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
            'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
            'post_process_flag',1,'err_cal_flag',1,'tol_adp',TOL_list(ii));
        end
        

        if flag_k3>0
            fprintf('--------- k = 3 --------------\n')
            main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
            'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
            'post_process_flag',1,'err_cal_flag',1,'tol_adp',TOL_list(ii));
        end
        

    end
 end

