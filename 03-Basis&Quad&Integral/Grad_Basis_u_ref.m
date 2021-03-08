function [Pr,Ps] = Grad_Basis_u_ref(a,b,i,j)
    % partial derivative of the basis_u

    h1 = GradJacobiP(a,0,0,i); 
    h2 = JacobiP(b,2*i+1,0,j);
    Pa = numeric_t('sqrt(2.0)')*h1.*h2.*(1-b).^i;

    g1 = JacobiP(a,0,0,i);
    g2 = JacobiP(b,2*i+1,0,j);
    g3 = GradJacobiP(b,2*i+1,0,j);

    Pb = numeric_t('sqrt(2.0)')*g1.*g3.*(1-b).^i ...
        + numeric_t('sqrt(2.0)')*g1.*g2.*(-i*(1-b).^(i-1));

    Pr = 2.0/(1-b)*Pa;
    Ps = (a+1)/(1-b) * Pa+Pb;
   
end