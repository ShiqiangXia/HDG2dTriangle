
%% define parameters

% problem
pb_type = 2011; % functional (u,g)
dom_ype = 'Rec';
Min_k = 2;
Max_k = 2;

% mesh
h0 = 0.1;
refine_flag = 1; % 1 adaptive; 0: uniform; -1: new uniform mesh each time
Niter_max = 1; % adaptive steps we do to locate trouble elements
Ncoarse_mesh = 3;% how many coarse mesh
N_outer_adap_steps = 5;

% other
tol_adp = 10e-14;
pp_flag = 1; % post_processing
err_cal_flag = 1;

% data
smooth_flag = 0; % 1: smooth data; 0: non-smooth
if smooth_flag == 0
    % non-smooth
    pri  = 1;
    adj  = 1;  
    fprintf('\nCase: u corner singularity, g= 1\n')
else
    pri  = 0; 
    % 2: x+y
    % 0: sin
    % 3: constant
    % 4: sin(2pi x)
    adj  = 0;
    fprintf('\nCase: u smooth \n')
    
end
%fprintf('Background mesh size h0 = %.2f\n\n',h0) 

%% call the test script

global class_t;
for order = Min_k:Max_k
    
    
    para = SetParameters_Ellipitc('order',order, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
        'pb_type',pb_type,'dom_type',dom_ype,'primal',pri,'adjoint',adj,...
        'post_process_flag',pp_flag,'err_cal_flag',err_cal_flag,'tol_adp',tol_adp);

    class_t=para.precision;

    %% Test 

    test_conv_adapt(para,Ncoarse_mesh, N_outer_adap_steps)
    
end