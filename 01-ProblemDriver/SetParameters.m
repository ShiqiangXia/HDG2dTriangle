function para = SetParameters()
    


    % define parameters
    para = Parameter();

    %% Problem parameters -------------------------------------------------
    pb_type = 1010;
    % pb_type: abcd
    % a: PDE-1 /Functional-2, b:source problem-0 or eigen problem-1,
    % c: PDE type (Poission-1), D:functional type (Vol-1, Bdry-2, Eig-3, Non-0)

    uexact = @(p)sin(pi*p(:,1)).*sin(pi*p(:,2));
    qexact = @(p)[-pi*cos(pi*p(:,1)).*sin(pi*p(:,2)),-pi*sin(pi*p(:,1)).*cos(pi*p(:,2))];

    %source_f = @(p)2*pi^2*sin(pi*p(:,1)).*(sin(pi*p(:,2)));
    source_f = @(p) 0*p(:,1);
    uD = @(p) 0*p(:,1)+5*p(:,1);
    uN = @(p) 0*p(:,1);

    vexact = @(p)sin(pi*p(:,1)).*sin(pi*p(:,2));
    pexact = @(p)[-pi*cos(pi*p(:,1)).*sin(pi*p(:,2)),-pi*sin(pi*p(:,1)).*cos(pi*p(:,2))];

    source_g = @(p)2*pi^2*sin(pi*p(:,1)).*(sin(pi*p(:,2)));
    vD = @(p) 0*p(:,1);
    vN = @(p) 0*p(:,1);

    para = para.SetPb(pb_type,...
        'uexact',uexact,'qexact',qexact,...
        'source_f',source_f,'uD',uD, 'uN',uN,...
        'vexact',vexact,'pexact',pexact,...
        'source_g',source_g,'vD',vD,'vN',vN);

    %% domain and mesh parameters----------------------------------------------

    structure_flag = 1;
    dom_type = 'Rec';
    h0 = 1;
    dirichlet_flag = ["bottom","top","left","right"];
    neuman_flag = [];
    %  dom_type and boundary names:   USE " " for boundary names, NOT ' '
    % 'Rec': ["bottom","top","left","right"]
    % 'L': ["bottom","left","right_low","right_high","top_low","top_high"]
    % 'Cir': ["Circle"]
    % 'CirHole': ["InnerCircle","OuterCircle"]

    x1 = -1;
    y1 = -1;
    x2 = 1;
    y2 = 1;
    tri_dir = 0;

    para = para.SetMesh(structure_flag,dom_type,h0,dirichlet_flag,neuman_flag,x1,y1,x2,y2,tri_dir);

    %% numerical method parameters---------------------------------------------

    order = 2;
    tau = numeric_t('1.0');
    post_process_flag = 0;

    para = para.SetNM(order,tau,post_process_flag);

    %% experiment parameters --------------------------------------------------
    precision = 'double';
    GQ_deg = 10;
    Niter = 1;
    refine_flag = 0;
    report_flag = 1;
    visualize_flag = 0;

    Neig = 5;
    Max_iter = 100;
    tol_eig = 1e-8;

    tol_adp = 1e-6;
    percent = 0.3;

    para = para.SetExp(precision,GQ_deg,Niter,refine_flag,report_flag,visualize_flag,...
        'Neig',Neig, 'Max_iter',Max_iter,'tol_eig',tol_eig,...
        'tol_adp',tol_adp,'percent',percent);

end
