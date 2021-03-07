function [a,b] = RStoAB(r,s)

% function [a,b] = RStoAB(r,s)
% Purpose : Transfer from (r,s) -> (a,b) coordinates in triangle

    Np = length(r); 
    a = zeros(Np,1,numeric_t);

    for n=1:Np
      if(s(n) ~= 1)
         a(n) = numeric_t('2.0')*(1+r(n))/(1-s(n))-numeric_t('1.0');
      else
          a(n) = numeric_t('-1.0'); 
      end
    end
    b = s;
    
end