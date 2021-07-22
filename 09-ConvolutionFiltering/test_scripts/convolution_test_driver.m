
%% define parameters

% problem
pb_type = 1011;
dom_ype = 'Rec';

% k
order = 1;
% mesh
h0 = 0.05;
refine_flag = 0; % 1 adaptive; 0: uniform
Niter_max = 1;

% other
tol_adp = 10e-14;
pp_flag = 1; % post_processing
err_cal_flag = 1;

% data
smooth_flag =1; % 1: smooth data; 0: non-smooth
if smooth_flag == 0
    pri  = 1;
    adj  = 1;  
    fprintf('\nCase: u corner singularity\n')
else
    pri  = 4; 
    % 2: x+y
    % 0: sin
    % 3: constant
    % 4: sin(2pi x)
    adj  = 0;
    fprintf('\nCase: u smooth \n')
end


%% set up parameters
global class_t;
para = SetParameters_Ellipitc('order',order, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',pri,'adjoint',adj,...
    'post_process_flag',pp_flag,'err_cal_flag',err_cal_flag,'tol_adp',tol_adp);

class_t=para.precision;

%% Test 


test_uh_pts(para)




