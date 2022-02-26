
%% define parameters

% problem
func_id = 2;

if func_id==1
    pb_type = 2011; % functional (u,g)
elseif func_id==2
    pb_type = 2012;
end
dom_ype = 'Rec';
Min_k = 1;
Max_k = 1;

% mesh
h0 = 0.1;
refine_flag = 1; % 1 adaptive; 0: uniform; -1: new uniform mesh each time

N_initial_iter_max = 1; % adaptive steps we do to locate trouble elements

Ncoarse_mesh = 3;% how many coarse mesh

N_outer_adap_steps = 8;

% other
tol_adp = 10e-14;

pp_flag = 1; % post_processing

err_cal_flag = 1;

% data
smooth_flag = 0; % 1: smooth data; 0: non-smooth
if smooth_flag == 0
    % non-smooth
    if func_id==1
        pri  = 1;
        adj  = 1;  
        fprintf('\nCase: u corner singularity, g= 1\n')
    elseif func_id==2
        pri  = 1;
        adj  = 2;
        fprintf('\n\nCase: u corner singularity and v discontinuous \n')
        
    end
    
    
else
    % 2: x+y
    % 0: sin
    % 3: constant
    % 4: sin(2pi x)
    if func_type==1
        pri  = 0; 
        adj  = 0;
        fprintf('\nCase: u smooth \n')
    elseif func_id==2
        pri  = 0;
        adj  = 4;
        fprintf('\n\nCase: u smooth and g smooth\n')
    
    end
    
end
%fprintf('Background mesh size h0 = %.2f\n\n',h0) 

%% call the test script

global class_t;
precision = 'double';
%precision = 'mp' ;
for order = Min_k:Max_k
    
    fprintf('k = %d\n', order);
    class_t = precision; %para.precision;
    para = SetParameters_Ellipitc('order',order, 'h0',h0, 'Niter',N_initial_iter_max, 'refine_flag', refine_flag,...
        'pb_type',pb_type,'dom_type',dom_ype,'primal',pri,'adjoint',adj,...
        'post_process_flag',pp_flag,'err_cal_flag',err_cal_flag,'tol_adp',tol_adp);

    

    %% Test 
    t_start = cputime;
    test_conv_adapt(para,Ncoarse_mesh, N_outer_adap_steps)
    %test_conv_adapt_with_extraction(para,Ncoarse_mesh, N_outer_adap_steps)
    % time 
    t_end = cputime - t_start;
    t_h = floor(t_end/3600);
    t_min = floor((t_end - t_h*3600)/60);
    t_sec = t_end - t_h*3600 - t_min*60 ;
    fprintf('Total time %d h %d min %.1f s\n',t_h,t_min,t_sec);
    
end