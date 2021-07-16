function [a,b] = RStoAB(r,s)

% function [a,b] = RStoAB(r,s)
% Purpose : Transfer from (r,s) -> (a,b) coordinates in triangle

% (r,s) is the coordinates on triagle (-1,-1), (1,-1) and (-1,1)
% (a,b) is the coordinates on square (-1,-1)x(1,1)

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