function ndof = GetDof(mymesh,k)
    
    % interior dof
    Nu = (k+2)*(k+1)/2;
    Nq = 2*Nu;
    Nuhat = k+1;
    
    N_interior = (Nu+Nq)*mymesh.num_elements;
    
    % skeleton dof
    % count interior faces
    N_skeleton = sum(mymesh.f_type==0) * Nuhat;
    
    ndof = N_interior+N_skeleton;
    
    
end