
h0 = 0.2;
refine_flag = 1;


fprintf('\n\nTest Eigenvalue problem\n\n')
fprintf('Initial mesh h0 = %f, tag_eig = 1;\n',h0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if refine_flag == 0
    fprintf('Uniform Refinement')
elseif refine_flag == 1
    fprintf('Adaptive Refinement\n')
end

fprintf('\n\nCase 1:  unit square\n')
Niter_max = 10;
pb_type = 2110;
dom_ype = 'Rec';


fprintf('--------- k = 1 --------------\n')
main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
    'post_process_flag',1,'err_cal_flag',1);

fprintf('--------- k = 2 --------------\n')
main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
    'post_process_flag',1,'err_cal_flag',1);

fprintf('--------- k = 3 --------------\n')
main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
    'post_process_flag',1,'err_cal_flag',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pb_type = 2110;
dom_ype = 'L';
Niter_max = 15;

fprintf('\n\nCase 2: L-shape domain \n')

fprintf('--------- k = 1 --------------\n')
main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
    'post_process_flag',1,'err_cal_flag',1);

fprintf('--------- k = 2 --------------\n')
main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
    'post_process_flag',1,'err_cal_flag',1);

fprintf('--------- k = 3 --------------\n')
main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',0,...
    'post_process_flag',1,'err_cal_flag',1);




