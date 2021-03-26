function [uh,qh,uh_hat] = HDG_SourcePbSolver_Elliptic(pb,mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,tau,source_f,uD,uN)
    
    % this source problem solver works for any source problem with dirichlet boundary
    % namely solve:
    % sum <qh*n+tau(uh-uhat),mu>_pK = 0
    %
    % Neuman boundary is not implemented yet. 
    %
    % for different type of PDE,the difference is to call different local solver 
    % matrices in Step 1. 
    % the rest assmeble global matrix, solve and recover part is the same!
    
    %% ---- step 1: Get local equations
    
    [List_LocSol, List_LocSol_f, List_Ns]...
        = HDG_GetLocalEquations_Elliptic(pb, mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,tau,source_f);
    
    %% ---- step 2: Assembel Global equations 
    
    [Global_M]= HDG_AssembleGlobalEquation_LHS_Elliptic(mymesh,...
        k,tau, List_LocSol,List_Ns);
    
    [Global_b]= HDG_AssembleGlobalEquation_RHS_Elliptic(mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,uD,uN, List_LocSol_f,List_Ns);
    
    
    %% ---- Step 3: Solve the linear system
    
    uh_hat = Global_M\Global_b;
    
    %% ---- Step 4: Recover the local variables. 
    
    [qh,uh] = HDG_RecoverLovalVariable_Elliptic(mymesh,...
        k,uh_hat,List_LocSol,List_LocSol_f);

    
    
end