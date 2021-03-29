function [lamh,wh_Neig,uh_Neig,ph_Neig,hat_var_Neig]...
        = HDG_EigPbSolver_Maxwell(pb,mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                    k,tau_t,tau_n,mu,epsilon,omg,Neig,Max_iter,Tol_eig)
    
      %% ---- step 1: Get local equations
      [List_LocSol, List_LocSol_f, List_Ns]...
        = HDG_Eig_GetLocalEquations_Maxwell(pb, mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,tau_t,tau_n,mu,epsilon,omg);
    
      %% ---- step 2: Assembel Global equations
      [Global_A]= HDG_AssembleGlobalEquation_LHS_Maxwell(mymesh,...
        k,tau_t,tau_n,epsilon, List_LocSol,List_Ns);
    
      %% ---- Step 3: Solve the eigenvalue problem iteratively
      
      [lamh,hat_var_Neig] = HDG_Eig_SolveEigenSystem_Maxwell(mymesh,k,...
          Global_A,List_Ns,List_LocSol_f,List_LocSol,...
          Neig,Max_iter,Tol_eig);
     
      
      %% ---- Step 4: Recover the local variables. 
      
      [wh_Neig,uh_Neig,ph_Neig] = HDG_Eig_RecoverLovalVariable_Maxwell(mymesh,...
        k,lamh,hat_var_Neig,List_LocSol,List_LocSol_f,Neig);

                
end