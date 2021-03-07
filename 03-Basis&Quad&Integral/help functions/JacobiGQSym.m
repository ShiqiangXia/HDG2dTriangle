function [x,w] = JacobiGQSym(alpha,beta,N);

% function [x,w] = JacobiGQ(alpha,beta,N)
% Purpose: Compute the N'th order Gauss quadrature points, x,
%          and weights, w, associated with the Jacobi
%          polynomial, of type (alpha,beta) > -1 ( <> -0.5).

if (N==0) x(1)=(alpha-beta)/(alpha+beta+2); w(1) = 2; return; end;
% Form symmetric matrix from recurrence.
J = sym(zeros(N+1));
diag2 = sym(zeros(N+1));

h1 = vpa(2*(0:N)+alpha+beta);
cc = vpa(2./(h1(1:N)+2).*sqrt((1:N).*((1:N)+alpha+beta).*...
    ((1:N)+alpha).*((1:N)+beta)./(h1(1:N)+1)./(h1(1:N)+3)));
diag2(1:N, 2:N+1) = diag(cc);
J = diag(vpa(-1/2*(alpha^2-beta^2)./(h1+2)./h1)) + ...
    diag2;
if (alpha+beta<10*eps) J(1,1)=0;end;
J = J + J';

%Compute quadrature by eigenvalue solve
J = vpa(J);
% [V,D] = eig(J); x = diag(D);
[V1,D] = eig(J); x1 = diag(D);
[x, Ind] = sort(x1,'ascend');
V = V1(:,Ind);

w = (V(1,:)').^2*2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
    gamma(beta+1)/gamma(alpha+beta+1);

end