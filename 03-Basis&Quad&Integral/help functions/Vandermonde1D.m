function [V1D] = Vandermonde1D(N,r)

% function [V1D] = Vandermonde1D(N,r)
% Purpose : Initialize the 1D Vandermonde Matrix, V_{ij} = phi_j(r_i);
% N: polynomial degree  
% r: points
% return matrix (dimension will be N+1 * N+1) 
% every row : same point different polynomials


V1D = zeros(length(r),N+1,numeric_t);
for j=1:N+1 % for each column
    V1D(:,j) = JacobiP(r(:), 0, 0, j-1);
    %V1D(:,j) = JacobiP(r(:), 0, 0, j-1)*sqrt(2/(2*j-1));
end
return