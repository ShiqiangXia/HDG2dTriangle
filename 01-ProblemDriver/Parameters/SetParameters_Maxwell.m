function para = SetParameters_Maxwell(varargin)
    
    [order,h0,Niter,refine_flag] = MyParaParse(varargin,'order','h0','Niter','refine_flag');

    % define parameters
    para = Parameter();

    %% Problem parameters -------------------------------------------------
    pb_type = 1150;
    % pb_type: abcd
    % a: PDE-1 /Functional-2, b:source problem-0 or eigen problem-1,
    % c: PDE type (Poission-1,Maxwell-5), D:functional type (Vol-1, Bdry-2, Non-0)

    %%%%%%%%%% smooth solution u %%%%%%%%
    mu = 1;
    epsilon = 1;
    omg = 0;
    
    % 1/mu - epsilon = 0 so that j = 0
    
    x = @(p) p(:,1);
    y = @(p) p(:,2);
    
    uexact_1 = @(p) sin(omg*y(p));
    uexact_2 = @(p) sin(omg*x(p));
    
    wexact = @(p) (omg/mu)*( cos(omg*x(p)) - cos(omg*y(p)) );
    
    pexact = @(p) 0*p(:,1);
    
    
    
    source_j_1 = @(p) 0*p(:,1);
    source_j_2 = @(p) 0*p(:,1);
   
    uxn_D = @(p,n) sin(omg*y(p))*n(2) -  sin(omg*x(p))*n(1);
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     mu = 1;
%     epsilon = 1;
%     omg = 1;
%     mypi = numeric_t('pi');
%     
%     x = @(p) p(:,1);
%     y = @(p) p(:,2);
%     
%     uexact_1 = @(p) (0*y(p));
%     uexact_2 = @(p) (0*x(p));
%     
%     wexact = @(p) 0*p(:,1);
%     
%     pexact = @(p) 0*p(:,1)+1;
%     source_j_1 = @(p) 0*p(:,1);
%     source_j_2 = @(p) 0*p(:,1);
%    
%     uxn_D = @(p,n) uexact_1(p)*n(2) -  uexact_2(p)*n(1);
%     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     mu = 1;
%     epsilon = 1;
%     omg = 1;
%     mypi = numeric_t('pi');
%     
%     x = @(p) p(:,1);
%     y = @(p) p(:,2);
%     
%     uexact_1 = @(p) (0*y(p)+1);
%     uexact_2 = @(p) (0*x(p)+1);
%     
%     wexact = @(p) 0*p(:,1);
%     
%     pexact = @(p) 0*p(:,1);
%     
%     source_j_1 = @(p) 0*p(:,1)-1;
%     source_j_2 = @(p) 0*p(:,1)-1;
%    
%     uxn_D = @(p,n) uexact_1(p)*n(2) -  uexact_2(p)*n(1);
%     
%     mu = 1;
%     epsilon = 1;
%     omg = 1;
%     mypi = numeric_t('pi');
%     
%     x = @(p) p(:,1);
%     y = @(p) p(:,2);
%     
%     uexact_1 = @(p) (y(p).^2);
%     uexact_2 = @(p) (x(p).^2);
%     
%     wexact = @(p) 2*x(p) - 2*y(p);
%     
%     pexact = @(p) 0*p(:,1);
%     
%     source_j_1 = @(p) -(y(p).^2+2);
%     source_j_2 = @(p) -(x(p).^2+2);
%    
%     uxn_D = @(p,n) uexact_1(p)*n(2) -  uexact_2(p)*n(1);
%     
%     
%     
   
    
    
    
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    para = para.SetPb(pb_type,...
        'wexact',wexact,'uexact_1',uexact_1,'uexact_2',uexact_2,'pexact',pexact,...
        'source_j_1',source_j_1,'source_j_2',source_j_2, 'uxn_D',uxn_D,...
        'mu',mu,'epsilon',epsilon,'omg',omg);

    %% domain and mesh parameters----------------------------------------------

    structure_flag = 0;
    
    dom_type = 'Rec';
    dirichlet_flag = ["bottom","top","left","right"];
    neuman_flag = [];
    x1 = 0;
    y1 = 0;
    x2 = 1;
    y2 = 1;
    tri_dir = 0;
    para = para.SetMesh(structure_flag,dom_type,h0,dirichlet_flag,neuman_flag,x1,y1,x2,y2,tri_dir);
%     

%     dom_type = 'L';
%     dirichlet_flag = ["bottom","top_high","right_low","left","right_high","top_low"];
%     neuman_flag = [];
%     tri_dir = 1;
%     para = para.SetMesh(structure_flag,dom_type,h0,dirichlet_flag,neuman_flag,tri_dir);
%     
    
    

%      
%     dom_type = 'CirHole'; 
%     dirichlet_flag = ["InnerCircle","OuterCircle"];
%     neuman_flag = [];
%     x1 = 0;
%     y1 = 0;
%     r1 = 1;
%     x2=0.2;
%     y2=0.2;
%     r2=0.5;
%     para = para.SetMesh(structure_flag,dom_type,h0,dirichlet_flag,neuman_flag,x1,y1,r1,x2,y2,r2);
%     
    
    %  dom_type and boundary names:   USE " " for boundary names, NOT ' '
    % 'Rec': ["bottom","top","left","right"]
    % 'L': ["bottom","left","right_low","right_high","top_low","top_high"]
    % 'Cir': ["Circle"]
    % 'CirHole': ["InnerCircle","OuterCircle"]


    

    %% numerical method parameters---------------------------------------------

    %order = 1;
    tau = numeric_t('1.0');
    post_process_flag = 0;

    para = para.SetNM(order,tau,post_process_flag);

    %% experiment parameters --------------------------------------------------
    precision = 'double';
    GQ_deg = 20;  % need more quads for corner singularity case
    
    %Niter = 3;
    %refine_flag = 1; 
    % 0: uniform refine,
    % -1: build new mesh based on h
    % 1: adaptive refine 'RGB', '2': 'RG' 3. 'NVB'
    
    err_cal_flag = 0; % 1: calculate L2 error of uh,qh
    
    report_flag = 1; 
    visualize_flag = 0;

    Neig = 3;
    Max_iter = 100;
    tol_eig = 1e-10;

    tol_adp = 1e-6;
    percent = 0.5;

    para = para.SetExp(precision,GQ_deg,Niter,refine_flag,err_cal_flag,report_flag,visualize_flag,...
        'Neig',Neig, 'Max_iter',Max_iter,'tol_eig',tol_eig,...
        'tol_adp',tol_adp,'percent',percent);

end
