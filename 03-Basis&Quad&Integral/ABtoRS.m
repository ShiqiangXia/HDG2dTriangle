function [r,s] = ABtoRS(a,b)
    r = (a+1)*(1-b)/2.0 - 1;
    s = b;
end
