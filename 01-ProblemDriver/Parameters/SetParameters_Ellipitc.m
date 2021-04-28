function para = SetParameters_Ellipitc(varargin)
    
    [order,h0,Niter,refine_flag,...
    pb_type,dom_type,...
    primal_data,adjoint_data,post_process_flag,err_cal_flag]...
    = MyParaParse(varargin,...
    'order','h0','Niter','refine_flag',...
    'pb_type','dom_type','primal','adjoint',...
    'post_process_flag','err_cal_flag');

    % define parameters
    para = Parameter();

    %% Problem parameters -------------------------------------------------
    %pb_type = 2111;
    % pb_type: abcd
    % a: PDE-1 /Functional-2, b:source problem-0 or eigen problem-1,
    % c: PDE type (Poission-1), D:functional type (Vol-1, Bdry-2, Non-0)


    % define the primal problem data
    if primal_data == 0 
        %%%%%%%%%% smooth solution u %%%%%%%%
        mypi = pi;
        uexact = @(p)sin(mypi*p(:,1)).*sin(mypi*p(:,2));
        qexact_1 = @(p)-mypi*cos(mypi*p(:,1)).*sin(mypi*p(:,2));
        qexact_2 = @(p)-mypi*sin(mypi*p(:,1)).*cos(mypi*p(:,2));
        source_f = @(p)2*mypi^2 * ( sin(mypi*p(:,1)).* sin(mypi*p(:,2)) );
        uD = uexact;
        uN = @(p) 0*p(:,1);
    
    elseif primal_data == 1
        %%%%%%%%%% cornor singularity solution u %%%%%%%%

        theta = @(p) atan(abs(p(:,2)./p(:,1))); % only consider 0<=theta<=pi/2
        radi  =  @(p) sqrt(p(:,1).^2+p(:,2).^2);

        uexact = @(p) radi(p).^(2/3) .* ( sin( 2/3*theta(p) ) );

        qexact_1 = @(p)-( 2/3 * p(:,1).* radi(p).^(-4/3).*(sin(2/3 *theta(p)))...
            - 2/3 * p(:,2).* radi(p).^(-4/3).*(cos(2/3 *theta(p))));
        qexact_2 = @(p) -(2/3 * p(:,2).* radi(p).^(-4/3).*(sin(2/3 *theta(p)))...
            + 2/3 * p(:,1).* radi(p).^(-4/3).*(cos(2/3 *theta(p))));

        source_f =  @(p) 0*p(:,1);

        uD = uexact;
        uN = @(p) 0*p(:,1);
    else
    %%%%%%%%%%%% Polynomial solution u %%%%%%%%%%%%%%%%%%%%
%     uexact = @(p) p(:,1)+p(:,2);
%     
%     qexact_1 = @(p) 0*p(:,1) - 1;
%     qexact_2 = @(p) 0*p(:,2) - 1;
%     source_f = @(p) 0*p(:,1);
%     uD = uexact;
%     uN = @(p) 0*p(:,1);
%     


    end 



    %--------------------------------------------------
    % define the adjoint problem data
    if adjoint_data == 0

    %%%%%%%%%% smooth solution v %%%%%%%%%%%%%%%%%%%%%%
        mypi = pi;
        vexact = @(p)sin(mypi*p(:,1)).*sin(mypi*p(:,2));
        pexact_1 = @(p)-mypi*cos(mypi*p(:,1)).*sin(mypi*p(:,2));
        pexact_2 = @(p)-mypi*sin(mypi*p(:,1)).*cos(mypi*p(:,2));

        source_g = @(p)2*mypi^2*sin(mypi*p(:,1)).*(sin(mypi*p(:,2)));
        vD = @(p) 0*p(:,1);
        vN = @(p) 0*p(:,1);
    
    elseif adjoint_data == 1
        %%%%%%%%%%% constant g = 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%

        vexact = @(p)0*p(:,1);
        pexact_1 = @(p)0*p(:,1);
        pexact_2 = @(p)0*p(:,1);

        source_g = @(p) 0*p(:,1)+1; % p(:,1).^4+2*p(:,2).^5; % 

        vD = @(p) 0*p(:,1);
        vN = @(p) 0*p(:,1);

    elseif adjoint_data == 2
   %%%%%%%%% v is discontinuous on the boundary %%%%%%%%%%%%%%%%%%%%%%
        aa = 0.2;
        bb = 0.8;
        condition_func =  @(p) (abs(p(:,2))<1e-10).*(p(:,1)>=aa).*(p(:,1)<=bb);
        vexact = @(p) condition_func(p).*( (p(:,1)-aa).*(p(:,1)-bb)+1);
        pexact_1 = @(p)0*p(:,1);
        pexact_2 = @(p)0*p(:,1);

        source_g = @(p)0*p(:,1);
        vD = vexact;
        vN = @(p) 0*p(:,1);
    elseif adjoint_data == 3
      %%%%%%%%%  v is continuous on the boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vexact = @(p)p(:,1).*p(:,2); % x*y
        pexact_1 = @(p)-p(:,2);
        pexact_2 = @(p)-p(:,x);

        source_g = @(p) 0*p(:,1);
        vD = @(p) p(:,1).*p(:,2);
        vN = @(p) 0*p(:,1);
    elseif adjoint_data == 4
        %%%%%%%%%  v smooth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mypi = pi;
        vexact = @(p)sin(mypi*p(:,1)).*exp(mypi*p(:,2)); % sin(pi x) *exp(pi y)
        pexact_1 = @(p)  - mypi * cos(mypi*p(:,1)).* exp(mypi*p(:,2));
        pexact_2  = @(p) - mypi * sin(mypi*p(:,1)).* exp(mypi*p(:,2));
        source_g = @(p) 0*p(:,1);
        vD = vexact;
        vN = @(p) 0*p(:,1);
        
    else

   
%     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     vexact = uexact;
%     pexact_1 = qexact_1;
%     
%     pexact_2 = qexact_2;
%     
%     source_g= source_f;
%     vD = vexact;
%     vN = @(p) 0*p(:,1);
%     
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
   
    
    
    tag_eig = 1;
    
    pb_text_info = WriteProbTextInfo(pb_type,primal_data,adjoint_data);

    para = para.SetPb(pb_type,...
        'uexact',uexact,'qexact_1',qexact_1,'qexact_2',qexact_2,...
        'source_f',source_f,'uD',uD, 'uN',uN,...
        'vexact',vexact,'pexact_1',pexact_1,'pexact_2',pexact_2,...
        'source_g',source_g,'vD',vD,'vN',vN,'tag_eig',tag_eig,...
        'pb_text_info',pb_text_info);

    %% domain and mesh parameters----------------------------------------------

    structure_flag = 1;
    
    %h0 = 0.5;
    if strcmp(dom_type,'Rec')
        %%%%% Rectangular %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dom_type = 'Rec';
        dirichlet_flag = ["bottom","top","left","right"];
        neuman_flag = [];
        x1 = 0;
        y1 = 0;
        x2 = 1;
        y2 = 1;
        tri_dir = 0;
        para = para.SetMesh(structure_flag,dom_type,h0,dirichlet_flag,neuman_flag,x1,y1,x2,y2,tri_dir);

    elseif strcmp(dom_type,'L')
        
        %%%%% L-shape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dom_type = 'L';
        dirichlet_flag = ["bottom","top_high","right_low","left","right_high","top_low"];
        neuman_flag = [];
        tri_dir = 0;
        para = para.SetMesh(structure_flag,dom_type,h0,dirichlet_flag,neuman_flag,tri_dir);

    elseif strcmp(dom_type,'CirHole')
             
        dom_type = 'CirHole'; 
        dirichlet_flag = ["InnerCircle","OuterCircle"];
        neuman_flag = [];
        x1 = 0;
        y1 = 0;
        r1 = 1;
        x2=0.2;
        y2=0.2;
        r2=0.5;
        para = para.SetMesh(structure_flag,dom_type,h0,dirichlet_flag,neuman_flag,x1,y1,r1,x2,y2,r2);
    else
        error('Wrong dom type! not implemented!')
    end
    
    
    
 
    %  dom_type and boundary names:   USE " " for boundary names, NOT ' '
    % 'Rec': ["bottom","top","left","right"]
    % 'L': ["bottom","left","right_low","right_high","top_low","top_high"]
    % 'Cir': ["Circle"]
    % 'CirHole': ["InnerCircle","OuterCircle"]


    

    %% numerical method parameters---------------------------------------------

    %order = 1;
    tau = numeric_t('1.0');
    % post_process_flag = 1
    para = para.SetNM(order,tau,post_process_flag);

    %% experiment parameters --------------------------------------------------

    precision = 'double';
    GQ_deg = 20;  % need more quads for corner singularity case
    
    %Niter = 3;
    %refine_flag = 1; 
    % 0: uniform refine,
    % -1: build new mesh based on h
    % 1: adaptive refine 'RGB', '2': 'RG' 3. 'NVB'
    
    %err_cal_flag = 1; % 1: calculate L2 error of uh,qh
    
    err_analysis_flag = 0; 
    % 1: compute the exact error terms like (q-qh,p-ph), need to know the exact solution
    
    report_flag = 1; 
    visualize_flag = 0;

    Neig = 3;
    Max_iter = 100;
    tol_eig = 1e-10;

    tol_adp = 1e-6;
    
    percent = 0.4;
    
    reduce_ratio=1e-10;

    para = para.SetExp(precision,GQ_deg,Niter,refine_flag,err_cal_flag,err_analysis_flag,report_flag,visualize_flag,...
        'Neig',Neig, 'Max_iter',Max_iter,'tol_eig',tol_eig,...
        'tol_adp',tol_adp,'percent',percent,'reduce_ratio',reduce_ratio);

end
