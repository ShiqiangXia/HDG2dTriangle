function [Auu,Auur,Auus,Auu3,Buuhat3] = HDG_LocalMatrix(k,GQ1DRef_pts,GQ1DRef_wts)
    % compute some basic local matrix on the reference triangle element
    Nu = (k+2)*(k+1)/2;
    
    
    % (basis_u1, basis_u2)
    Auu = eye(Nu,Nu,numeric_t);
    
    % Auur =  (basis_u_i, d(basis_u_j)/dr)
    % Auus = (basis_u_i, d(basis_u_j)/ds)
    
    % size Nu x Nu
    [Auur,Auus] = Volume_Int_u_du(k,GQ1DRef_pts,GQ1DRef_wts);
    
    % Auu3 = <basis_u1, basis_u2>_{F_i} i = 1,2,3
    % size Nu x Nu x 3
    Auu3 = Boundary_Int_u_u(k,GQ1DRef_pts,GQ1DRef_wts);
    
    %Buuhat3 =  <basis_u_i, basis_uhat_j>_{F_i} i = 1,2,3
    % Size Nu x Nuhat x 3
    Buuhat3 = Boundary_Int_u_uhat(k,GQ1DRef_pts,GQ1DRef_wts);
    
    
    
    
    
    
    
    
    
end