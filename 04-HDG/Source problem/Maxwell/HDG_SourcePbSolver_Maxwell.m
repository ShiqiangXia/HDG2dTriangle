function [wh,uh,ph,hat_var] = HDG_SourcePbSolver_Maxwell(pb,mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,tau_t,tau_n,source_j_1,source_j_2,uhat_t_D,mu,epsilon,omg)
    
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
        = HDG_GetLocalEquations_Maxwell(pb, mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,tau_t,tau_n,source_j_1,source_j_2,mu,epsilon,omg);
    
    %% ---- step 2: Assembel Global equations 
    
    [Global_M]= HDG_AssembleGlobalEquation_LHS_Maxwell(mymesh,...
        k,tau_t,tau_n,epsilon, List_LocSol,List_Ns);
    
    [Global_b]= HDG_AssembleGlobalEquation_RHS_Maxwell(mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,uhat_t_D, List_LocSol_f,List_Ns);
    
    
    %% ---- Step 3: Solve the linear system
    
    hat_var = Global_M\Global_b;
    
    %% ---- Step 4: Recover the local variables. 
    
    [wh,uh,ph] = HDG_RecoverLovalVariable_Maxwell(mymesh,...
        k,hat_var,List_LocSol,List_LocSol_f);

    
    
end