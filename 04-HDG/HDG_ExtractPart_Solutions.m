function [uh,qh,uhat] = HDG_ExtractPart_Solutions(k_in,k_out,uhstar,qhstar,uhatstar)
    % Extract HDG solution coefficients of lower order polynomials
    % input: HDG solution of degree k_in
    % output: extract the coefficients of degree k_out
    Nk_in =(k_in+1)*(k_in+2)/2;
    r_ind = Get_Hierarchical_indx(k_in,k_out);
    
    uh = uhstar(r_ind,:);

    qh = [qhstar(r_ind,:);qhstar(r_ind+Nk_in,:)];
    
    temp = reshape(uhatstar,k_in+1,[]);
    temp = temp(1:k_out+1,:);
    
    uhat = reshape(temp,[],1);
    
end