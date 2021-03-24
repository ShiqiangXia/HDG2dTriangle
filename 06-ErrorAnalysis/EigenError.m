function rlt = EigenError(pde_type,lamh,dom_type,geo_para)
    % step 1 get exact eigenvalues
    % step 2 compute error
    
    Neig = size(lamh,2);
    rlt = zeros(1,Neig,numeric_t);
    if strcmp(pde_type,'1')
        % poisson problem
        
        % rectangular domain, eigenvalue is pi^2 *( (m/Lx)^2+(n/Ly)^2 )
        if strcmp(dom_type,'Rec')
            x1 = geo_para{1};
            y1 = geo_para{2};
            x2 = geo_para{3};
            y2 = geo_para{4};
            Lx = abs(x1-x2);
            Ly = abs(y1-y2);
            N_eig_fix = 100;
            l = zeros(1,N_eig_fix*N_eig_fix,numeric_t);
            for m=1:N_eig_fix
                for n = 1:N_eig_fix
                    l(1,(m-1)*N_eig_fix+n) = (m/Lx)^2+(n/Ly)^2;
                end
            end
            l = sort(l);
            lam_exact = l*numeric_t('pi^2');
            % the first 885573 eigenvalues are corrrect!
            
            
            rlt(1,:) = abs(lam_exact(1:Neig) - lamh );
            
        elseif strcmp(dom_type,'L')
            
            if (Neig>3)
                Neig = 3;
                frpintf("L-shaped domain only implements the first 3 eigenvalues")
            end
            
            lam_exact = numeric_t('[9.6397238440219, 15.19725192576365 ,2*pi^2]');
            
            rlt(1,:) = abs(lam_exact(1:Neig) - lamh' );
            
        else
            error('Domain type %s has not been implemented yet',dom_type);
            
        end
    else
        error('PDE type %s has not been implemented yet',pde_type );
    end
    
   
    
end