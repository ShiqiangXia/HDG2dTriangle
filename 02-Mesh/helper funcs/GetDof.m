function ndof = GetDof(mymesh,k)
    
    % interior dof
    Nu = (k+2)*(k+1)/2;
    Nq = 2*Nu;
    
    N_interior = (Nu+Nq)*mymesh.num_elements;
    
    % skeleton dof
    % count interior faces
    mymesh.f_type
    
    
end