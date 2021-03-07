function [x] = JacobiGLSym(alpha,beta,N)

% function [x] = JacobiGLSym(alpha,beta,N)
% Purpose: Compute the N'th order Gauss Lobatto quadrature 
%          points, x, associated with the Jacobi polynomial,
%          of type (alpha,beta) > -1 ( <> -0.5). 

% alpha = sym(alpha);
% beta  = sym(beta);

x = sym(zeros(N+1,1));
if (N==0) x(1)=0.0; return; end;
if (N==1) x(1)=-1.0; x(2)=1.0; return; end;

[xint,w] = JacobiGQSym(alpha+1,beta+1,N-2);
x(1) = -1;  
x(2:N,1) = xint; 
x(N+1)=1;
% x = [-1, xint', 1]';
end



