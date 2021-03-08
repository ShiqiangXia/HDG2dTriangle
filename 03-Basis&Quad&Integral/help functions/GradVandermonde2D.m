function [V2Dr,V2Ds] = GradVandermonde2D(N,a,b)

% function [V2Dr,V2Ds] = GradVandermonde2D(N,a,b)
% (a,b) is the coordinates on the ref square
% Purpose : Initialize the gradient of the modal basis (i,j) 
% at points (a,b) 	

% every row : same point different polynomials

V2Dr = zeros(length(a),(N+1)*(N+2)/2); V2Ds = zeros(length(a),(N+1)*(N+2)/2);

% % find tensor-product coordinates
% [a,b] = RStoAB(a,b);

% Initialize matrices
sk = 1;
for i=0:N
  for j=0:N-i
    [V2Dr(:,sk),V2Ds(:,sk)] = Grad_Basis_u_ref(a,b,i,j);
    sk = sk+1;
  end
end

return;