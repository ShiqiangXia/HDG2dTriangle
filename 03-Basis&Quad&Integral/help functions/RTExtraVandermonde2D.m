function [V2D_1,V2D_2] = RTExtraVandermonde2D(N, a, b)

dimRT_extra = N+1;

V2D_1 = zeros(length(a),dimRT_extra,numeric_t);
V2D_2 = zeros(length(a),dimRT_extra,numeric_t);

for i = 0:N
    V2D_1(:,1+i) = Basis_u_ref(a,b,i+1,N-i);
    V2D_2(:,1+i) = Basis_u_ref(a,b,i,N+1-i);
end





return;