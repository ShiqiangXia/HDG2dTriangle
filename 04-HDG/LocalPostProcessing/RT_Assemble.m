function RT_mat = RT_Assemble()
    % assemble the post-processing matrix for each element
    %dimRT = (k+3)(k+1) 
    %      = 2*(k*(k+1)/2) + 3*(k+1)     test: P_{k-1} +sum P_{k}(F)
    %      = 2*(k+1)*(k+2)/2 + (k+1)     trial:P_{k}+(x,y)*P_k_homo
    % (qh*,tau)_K = Jk (qh*,tau)_{K_ref}
    
    %<qh*n,mu>_F = 0.5*|F| <qh*n,mu>_[-1,1]
    
    % ii*(k+1)+1: (ii+1)*(k+1), 1:(k+1)*(k+2)/2   <n1*pk,mu>
    
    % ii*(k+1)+1: (ii+1)*(k+1), (k+1)*(k+2)/2+1: (k+1)*(k+2) <n2*pk,mu>
    
    % ii*(k+1)+1: (ii+1)*(k+1), (k+1)*(k+2)+1:end   <n1*RT_ext_1,mu> + <n2*RT_ext_2,mu>
    
end