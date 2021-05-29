function type = CheckType(ii,l)
    
    d1 = l(ii-1); 
    
    d2 = l(ii);
    
    d3 = l(ii+1);
    
    test_1 = d1>d2;change_1 = (d1-d2)/d1 >0.05;
    test_2 = d2>d3;change_2 = (d2-d3)/d2 >0.05;
    
    if test_1 && test_2 && change_1 && change_2
        type = 1;
    elseif test_1 && change_1 &&  (~test_2 ||~change_2)
        type = 2;
    elseif (~test_1 || ~change_1) && test_2 && change_2
        type = 3;
    else
        type = 4;
    end
    
end