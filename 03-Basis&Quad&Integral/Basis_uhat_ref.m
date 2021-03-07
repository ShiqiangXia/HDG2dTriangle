function P = Basis_uhat_ref(xs,N)
    % Basis function for uhat  on reference element [-1,1]
    % normalized legendre polynomial
    % N>=1---> degree N-1 legendre polynomial
    if N<1
        error('Basis index N should >= 1')
    end
    P = JacobiP(xs,0,0,N-1); 
    
end
