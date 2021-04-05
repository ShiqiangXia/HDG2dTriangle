function [V2D_1,V2D_2] = RTVandermonde2D(N, a, b)

% function [V2D] = Vandermonde2D(N, a, b);
% (a,b) is the coordinates on the ref square
% Purpose : Initialize the 2D Vandermonde Matrix 
% for the Ravart-Thomas space of N, 
% Size #points x #basis
% RT basis (pk,0), (0,pk),  (p_{i+1,N-i}, p_{i,N+1-i}) i=0,1,...,N

% every row : same point different polynomials

dimRT =(N+3)*(N+1);
dimPk = (N+1)*(N+2)/2;
V2D_1 = zeros(length(a),dimRT,numeric_t);
V2D_2 = zeros(length(a),dimRT,numeric_t);

V2D = zeros(length(a),dimPk,numeric_t);

% % Transfer to (a,b) coordinates
% [a, b] = RStoAB(r, s);

% build the Vandermonde matrix
count_basis = 1;

for i=0:N
  for j=0:N - i
    V2D(:,count_basis) = Basis_u_ref(a,b,i,j);
    count_basis = count_basis+1;
  end
end

V2D_1(:,1:dimPk) = V2D;
V2D_2(:,dimPk+1:2*dimPk) = V2D;

for i = 0:N
    V2D_1(:,2*dimPk+1+i) = Basis_u_ref(a,b,i+1,N-i);
    V2D_2(:,2*dimPk+1+i) = Basis_u_ref(a,b,i,N+1-i);
end





return;