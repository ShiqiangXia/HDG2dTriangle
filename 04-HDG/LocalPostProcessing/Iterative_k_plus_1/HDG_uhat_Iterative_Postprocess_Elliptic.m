function [uh_out, qh_out, uhat_out]=HDG_uhat_Iterative_Postprocess_Elliptic(pb,mymesh,...
        GQ1DRef_pts,GQ1DRef_wts,...
        k_out,k_in,tau,uhat,source_f,uD,uN)
    
    % Iterative postprocessing
    % k_in= k; k_out = k+1
    % input: uhat_k from HDG k 
    % step 1: uhat_k+1 = Proj(uhat_k) (project uhat_k to P_{k+1} )
    % step 2: Solve HDG k+1 Local solver A_{k+1} x = uhat_k+1
    % step 3: get new uhat_k+1 using the relation of uhat and uh,qh
    %        formula uhat = tau+/T uh+  + tau-/T + uh-  + 1/T [[qh]]
    %         here T = tau+ + tau-,  [[qh]] = qh+ * n+  + qh- * n-
    % step 4: GO TO step 2 and iterate
    
    Max_iter = 0;
    
    %% ---- step 1: Project uhat_k to P_{k+1} space ----------------------
    % 1. set uhat to the right position
    % 2. for bdry, project uD to k+1
    [uhat0, uhatD] = ProjectUhat(mymesh, k_out, k_in, uhat, uD, uN, GQ1DRef_pts,GQ1DRef_wts);
    
    %% ---- step 2.1: Get local equations for HDG k_out -------------------
    
    [List_LocSol, List_LocSol_f, List_Ns]...
        = HDG_GetLocalEquations_Elliptic(pb, mymesh,GQ1DRef_pts,GQ1DRef_wts,...
        k_out,tau,source_f);
    
    %% ---- Step 2.2: Recover the local variables -------------------------
    [qh0,uh0] = HDG_RecoverLovalVariable_Elliptic(mymesh,...
        k_out,uhat0,List_LocSol,List_LocSol_f);
    
    qh_out = qh0 ;
    uh_out = uh0 ;
    uhat_out = uhat0;
    
     %% ---- Step 3: Iterative process -------------------------------------
    if Max_iter > 0
        for ii = 1:Max_iter

            uhat_out = BuildUhatFromQhUh(mymesh, k_out, qh_out, uh_out, tau, uhatD, List_Ns) ;
            [qh_out,uh_out] = HDG_RecoverLovalVariable_Elliptic(mymesh,...
            k_out,uhat_out,List_LocSol,List_LocSol_f);

        end
    end
    
    
    
    
end