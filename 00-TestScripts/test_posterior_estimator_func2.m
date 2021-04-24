
h0 = 0.5;
refine_flag = 1;
pb_type = 2012;

fprintf('\n\nTest Linear functinals J(u) = <q*n,psi>\n\n')
fprintf('Initial mesh h0 = %f\n',h0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if refine_flag == 0
    fprintf('Uniform Refinement')
elseif refine_flag == 1
    fprintf('Adaptive Refinement\n')
end

fprintf('\n\nCase 1: u smooth and v smooth on unit square\n')
Niter_max = 15;

dom_ype = 'Rec';


fprintf('--------- k = 1 --------------\n')
main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
    'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',4,...
    'post_process_flag',1,'err_cal_flag',1);

% fprintf('--------- k = 2 --------------\n')
% main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',4,...
%     'post_process_flag',1,'err_cal_flag',1);
% 
% fprintf('--------- k = 3 --------------\n')
% main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',4,...
%     'post_process_flag',1,'err_cal_flag',1);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% dom_ype = 'Rec';
% Niter_max = 15;
% 
% fprintf('\n\nCase 2: u smooth and v discontinuous \n')
% 
% fprintf('--------- k = 1 --------------\n')
% main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',2,...
%     'post_process_flag',1,'err_cal_flag',1);
% 
% fprintf('--------- k = 2 --------------\n')
% main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',2,...
%     'post_process_flag',1,'err_cal_flag',1);
% 
% fprintf('--------- k = 3 --------------\n')
% main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',0,'adjoint',2,...
%     'post_process_flag',1,'err_cal_flag',1);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% dom_ype = 'Rec';
% 
% fprintf('\n\nCase 3: u corner singularity and v discontinuous \n')
% 
% fprintf('--------- k = 1 --------------\n')
% main('elliptic','order',1, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',1,'adjoint',2,...
%     'post_process_flag',1,'err_cal_flag',1);
% 
% fprintf('--------- k = 2 --------------\n')
% main('elliptic','order',2, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',1,'adjoint',2,...
%     'post_process_flag',1,'err_cal_flag',1);
% 
% fprintf('--------- k = 3 --------------\n')
% main('elliptic','order',3, 'h0',h0, 'Niter',Niter_max, 'refine_flag', refine_flag,...
%     'pb_type',pb_type,'dom_type',dom_ype,'primal',1,'adjoint',2,...
%     'post_process_flag',1,'err_cal_flag',1);
% 

