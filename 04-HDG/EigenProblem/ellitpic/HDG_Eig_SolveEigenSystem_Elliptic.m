function [lamh,uhat_Neig] = HDG_Eig_SolveEigenSystem_Elliptic(mymesh,k,...
          Global_A,List_Ns,List_LocSol_f,List_LocSol,...
          Neig,Max_iter,Tol_eig)
      
      
      Nuhat = k+1;
      num_faces = mymesh.num_faces;
      
      lamh = zeros(1,Neig,numeric_t);
      uhat_Neig = zeros(num_faces*Nuhat,Neig,numeric_t);
      
      lam = 0;
      
      [M0,~] = HDG_Eig_AssembleMlambda_Elliptic(lam,mymesh,k,List_Ns,List_LocSol,List_LocSol_f); 
      
      [V0,D0] = eigs(Global_A,M0,Neig ,'sm');
      
      for ii = 1:Neig
          eta0 = V0(:,ii);
          lam0 = D0(ii,ii);
          max_iter_flag = 1;
          
          for jj = 1:Max_iter
            
            [Mlam,Nlam] = HDG_Eig_AssembleMlambda_Elliptic(lam0,mymesh,k,List_Ns,List_LocSol, List_LocSol_f); 

            eta_hat = (Global_A-lam0*Mlam)\(Nlam*eta0);

            delta = numeric_t('1.0')/(eta_hat'*eta0);

            lam = lam0 + delta;
            eta = delta*eta_hat;

            if abs(delta) < Tol_eig
                max_iter_flag = 0;
                break;
            else
                lam0 = lam;
                eta0 = eta;
            end
            
          end
          
        if(max_iter_flag==1)
            fprintf("%d th eigenvalue, Reached the max iteration: %d\n",ii, Max_iter)
        end

        lamh(1,ii) =lam;
        uhat_Neig(:,ii) = eta;
        
      end
      
end