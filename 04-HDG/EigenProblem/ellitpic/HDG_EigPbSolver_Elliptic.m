function [lamh,uh_Neig,qh_Neig,uhat_Neig]...
        = HDG_EigPbSolver_Elliptic(pb,mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                    k,tau,Neig,Max_iter,Tol_eig)
    
      %% ---- step 1: Get local equations
      [List_LocSol, List_LocSol_f, List_Ns]...
        = HDG_Eig_GetLocalEquations_Elliptic(pb, mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k,tau);
    
      %% ---- step 2: Assembel Global equations
      [Global_A]= HDG_AssembleGlobalEquation_LHS_Elliptic(mymesh,...
        k,tau, List_LocSol,List_Ns);
    
      %% ---- Step 3: Solve the eigenvalue problem iteratively
      
      [lamh,uhat_Neig] = HDG_Eig_SolveEigenSystem_Elliptic(mymesh,k,...
          Global_A,List_Ns,List_LocSol_f,List_LocSol,...
          Neig,Max_iter,Tol_eig);
     
      
      %% ---- Step 4: Recover the local variables. 
      
      [qh_Neig,uh_Neig] = HDG_Eig_RecoverLovalVariable_Elliptic(mymesh,...
        k,lamh,uhat_Neig,List_LocSol,List_LocSol_f,Neig);

                
end