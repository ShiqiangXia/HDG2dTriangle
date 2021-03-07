function [V2D] = Vandermonde2D(N, r, s)

% function [V2D] = Vandermonde2D(N, r, s);
% Purpose : Initialize the 2D Vandermonde Matrix, 
% Size #points x #basis
% V_{ij} = phi_(i,j)(r_t, s_t);

% every row : same point different polynomials

V2D = zeros(length(r),(N+1)*(N+2)/2,numeric_t);

% Transfer to (a,b) coordinates
[a, b] = RStoAB(r, s);

% build the Vandermonde matrix
count_basis = 1;

for i=0:N
  for j=0:N - i
    V2D(:,count_basis) = Basis_u_ref(a,b,i,j);
    count_basis = count_basis+1;
  end
end
return;