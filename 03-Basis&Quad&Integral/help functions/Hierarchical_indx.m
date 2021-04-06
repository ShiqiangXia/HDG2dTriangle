function r_ind = Hierarchical_indx(k)
    % get the hierarchical index of P_{k-1} in P_{k}
    % notice that P_{k} is order as 
    % i = 0,1,...k
    % j = 0,1,..k-i
    %  ->  P_{i,j}
    
    
    
    Nk_1 =  (k+1)*k/2;
    r_ind = zeros(Nk_1,1,numeric_t);
    
    ct = 1;
    sk=1;
    for ii = 0:k
        for jj = 0:k-ii
            if(ii+jj<k)
                r_ind(ct,1) = sk;
                ct = ct+1;
            end
            sk = sk+1;
        end
    end
    
end