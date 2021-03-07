function P=Basis_u_ref(a,b,i,j)
    
    % basis function for scalar variable
    % function [P] = Basis_u_ref(a,b,i,j);
    % 
    % Purpose : Evaluate 2D orthonormal polynomial
    %           on triangle of order (i,j).
    % Dubiner basis!
    % Here (a,b) are coordinates on reference square [-1,1]^2
    %

    h1 = JacobiP(a,0,0,i); 
    h2 = JacobiP(b,2*i+1,0,j);
    P = numeric_t('sqrt(2.0)')*h1.*h2.*(1-b).^i; 
    
    
end