% Test reduce ratio results

% reduce_ratio = 1e-4;
% adptive mark strategy: esti_K > percent * max(esti_K);
% percent = 0.5;
% intial mesh h = 0.2;

Niter_max = 100;
h0 = 1;
fprintf('Test reduce_ratio = 1e-6; percent = 0.5; intial mesh h = %.2f\n',h0)

fprintf('\n\nTest Linear functinals J(u) = (u,g)\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\nCase 1: u smooth and g smooth on unit square\n')

% pb_type = 2011;
% dom_ype = 'Rec';
% 
% fprintf('--------- k = 1 --------------\n')
% main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', 1,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
%     'post_process_flag',0,'err_cal_flag',1);
% 
% fprintf('--------- k = 2 --------------\n')
% main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', 1,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
%     'post_process_flag',0,'err_cal_flag',1);
% 
% fprintf('--------- k = 3 --------------\n')
% main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', 1,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
%     'post_process_flag',0,'err_cal_flag',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pb_type = 2011;
dom_ype = 'Rec';

fprintf('\n\nCase 2: u corner singularity and g=1 smooth\n')

fprintf('--------- k = 1 --------------\n')
main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', 1,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',1,'adjoint',1,...
    'post_process_flag',0,'err_cal_flag',1);

fprintf('--------- k = 2 --------------\n')
main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', 1,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',1,'adjoint',1,...
    'post_process_flag',0,'err_cal_flag',1);

fprintf('--------- k = 3 --------------\n')
main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', 1,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',1,'adjoint',1,...
    'post_process_flag',0,'err_cal_flag',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n\nTest Eigenvalue problems, target 1st eigenvalue\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\nCase 1: unit square domain\n')
% 
% pb_type = 2111;
% dom_ype = 'Rec';
% 
% fprintf('--------- k = 1 --------------\n')
% main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', 1,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
%     'post_process_flag',0,'err_cal_flag',0);
% 
% fprintf('--------- k = 2 --------------\n')
% main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', 1,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
%     'post_process_flag',0,'err_cal_flag',0);
% 
% fprintf('--------- k = 3 --------------\n')
% main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', 1,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
%     'post_process_flag',0,'err_cal_flag',0);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\nCase 2: L-shape domain\n')

pb_type = 2111;
dom_ype = 'L';
fprintf('--------- k = 1 --------------\n')
main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', 1,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
    'post_process_flag',0,'err_cal_flag',0);

fprintf('--------- k = 2 --------------\n')
main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', 1,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
    'post_process_flag',0,'err_cal_flag',0);

fprintf('--------- k = 3 --------------\n')
main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', 1,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
    'post_process_flag',0,'err_cal_flag',0);


