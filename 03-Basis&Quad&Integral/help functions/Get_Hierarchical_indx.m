function r_ind = Get_Hierarchical_indx(k_in,k_out)
    % get the hierarchical index of P_{k_out} in P_{k_in}
    % notice that P_{k} is order as 
    % i = 0,1,...k
    % j = 0,1,..k-i
    %  ->  P_{i,j}
    
    if k_out>k_in
        error('Wrong inputs! we need k_out <= k_in')
    end
    
    Nk_out =  (k_out+2)*(k_out+1)/2;
    
    r_ind = zeros(Nk_out,1,numeric_t);
    
    ct = 1;
    sk=1;
    
    for ii = 0:k_in
        for jj = 0:k_in-ii
            if(ii+jj<=k_out)
                r_ind(ct,1) = sk;
                ct = ct+1;
            end
            sk = sk+1;
        end
    end
    
end