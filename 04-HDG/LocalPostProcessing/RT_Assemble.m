function RT_mat = RT_Assemble(k,Jk,vertice_list, A_qtau,Buuhat3,RTuuhat3)
    % assemble the post-processing matrix for each element
    %dimRT = (k+3)(k+1) 
    %      = 2*(k*(k+1)/2) + 3*(k+1)     test: P_{k-1} +sum P_{k}(F)
    %      = 2*(k+1)*(k+2)/2 + (k+1)     trial:P_{k}+(x,y)*P_k_homo
    
    
    dimRT = (k+3)*(k+1);
    Nk_1 = (k*(k+1))/2;
    Nk = (k+1)*(k+2)/2;
    Nmu = k+1;
    
    RT_mat = zeros(dimRT,dimRT,numeric_t);
    % (qh*,tau)_K = Jk (qh*,tau)_{K_ref}
    RT_mat(1:Nk_1,1:Nk) = Jk*A_qtau;
    RT_mat(Nk_1+1:2*Nk_1,Nk+1:2*Nk) = Jk*A_qtau;
    
    %<qh*n,mu>_F = 0.5*|F| <qh*n,mu>_[-1,1]
    
    % ii*(k+1)+1: (ii+1)*(k+1), 1:(k+1)*(k+2)/2   <n1*pk,mu>
    
    % ii*(k+1)+1: (ii+1)*(k+1), (k+1)*(k+2)/2+1: (k+1)*(k+2) <n2*pk,mu>
    
    % ii*(k+1)+1: (ii+1)*(k+1), (k+1)*(k+2)+1:end   <n1*RT_ext_1,mu> + <n2*RT_ext_2,mu>
    
    [e_list,n1,n2,n3] = GetTriFaceInfo(vertice_list);
    for ii = 1:3
        
        if ii == 1
            norm_vec = n1;
        elseif ii == 2
            norm_vec = n2;
        else
            norm_vec = n3;
        end
        
        start_id = 2*Nk_1 + (ii-1)*Nmu+1;
        end_id = start_id-1+Nmu;
        
        RT_mat(start_id:end_id,1:Nk) = 0.5*e_list(ii)*norm_vec(1)*(Buuhat3(:,:,ii))';
        RT_mat(start_id:end_id,Nk+1:2*Nk) = 0.5*e_list(ii)*norm_vec(2)*(Buuhat3(:,:,ii))';
        
        RT_mat(start_id:end_id,2*Nk+1:end) = 0.5*e_list(ii)...
                                             * (norm_vec(1)*RTuuhat3(:,:,ii,1)'...
                                             +norm_vec(2)*RTuuhat3(:,:,ii,2)' );
        
    end
    
end