function type = CheckType(ii,l)
    
    d1 = l(ii-1) - l(ii);
    
    d2 = l(ii) - l(ii+1);
    
    if d1>0 && d2>0
        type = 1;
    elseif d1>0 && d2<=0
        type = 2;
    elseif d1<=0 && d2>0
        type = 3;
    else
        type = 4;
    end
    
end