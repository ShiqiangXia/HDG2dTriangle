% script to test solve PDE with HDG method

refine_flag = 0; % 1 adaptive; 0: uniform

Niter_max = 3;

%TOL_list = [1e-14];
%N_TOL = size(TOL_list,2);
tol_adp = 10e-14;

smooth_flag =1; % 1: smooth data; 0: non-smooth
pp_flag = 1; % post_processing
err_cal_flag = 1; 

% 1: run this case; 0: not run

flag_k1    = 1; 
flag_k2    = 1;
flag_k3    = 1;

if refine_flag == 0
    fprintf('Uniform Refinement\n')
elseif refine_flag == 1
    fprintf('Adaptive Refinement\n')
end

h0 = 0.2;

fprintf('--------------------------------------------------------------\n')
fprintf('--------------------------------------------------------------\n')
fprintf('Initial mesh h0 = %f\n',h0);
fprintf('Solve PDE with HDG method \n')


pb_type = 1011;
dom_ype = 'Rec';

if smooth_flag == 0
    pri  = 1;
    adj  = 1;  
    fprintf('\nCase: u corner singularity\n')
else
    pri  = 0; % 2
    adj  = 0;
    fprintf('\nCase: u smooth \n')
end


if flag_k1>0
    fprintf('--------- k = 1 --------------\n')
    main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',pri,'adjoint',adj,...
    'post_process_flag',pp_flag,'err_cal_flag',err_cal_flag,'tol_adp',tol_adp);
end


if flag_k2>0
    fprintf('--------- k = 2 --------------\n')
    main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',pri,'adjoint',adj,...
    'post_process_flag',pp_flag,'err_cal_flag',err_cal_flag,'tol_adp',tol_adp);
end

if flag_k3>0
    fprintf('--------- k = 3 --------------\n')
    main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',pri,'adjoint',adj,...
    'post_process_flag',pp_flag,'err_cal_flag',err_cal_flag,'tol_adp',tol_adp);

end
