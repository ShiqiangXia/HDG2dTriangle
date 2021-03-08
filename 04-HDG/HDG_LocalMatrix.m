function [Aqq,Auur,Auus,Auu3,Buuhat3] = HDG_LocalMatrix(k,GQ1DRef_pts,GQ1DRef_wts)
    % compute some basic local matrix on the reference triangle element
    Nu = (k+2)*(k+1)/2;
    Nq = 2*Nu;
    
    % (basis_q1, basis_q2)
    Aqq = eye(Nq,Nq,numeric_t);
    
    % (basis_u, d(basis_u)/dr)
    % (basis_u, d(basis_u)/ds)
    % size Nu x Nu
    [Auur,Auus] = Volume_Int_u_du(k,GQ1DRef_pts,GQ1DRef_wts);
    
    % <basis_u1, basis_u2>_{F_i} i = 1,2,3
    % size Nu x Nu x 3
    Auu3 = Boundary_Int_u_u(k,GQ1DRef_pts,GQ1DRef_wts);
    
    % <basis_u, basis_uhat>_{F_i} i = 1,2,3
    % Size Nu x Nuhat x 3
    Buuhat3 = Boundary_Int_u_uhat(k,GQ1DRef_pts,GQ1DRef_wts);
    
    
    
    
    
    
    
    
    
end