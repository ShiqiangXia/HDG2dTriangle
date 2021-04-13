function [uh,qh,uhat] = HDG_ExtractPart_Solutions(k_in,k_out,uhstar,qhstar,uhatstar)
    % Extract HDG solution coefficients of lower order polynomials
    % input: HDG solution of degree k_in
    % output: extract the coefficients of degree k_out
    
    r_ind = Get_Hierarchical_indx(k_in,k_out);
    
    uh = uhstar(r_ind,:);

    qh = [qhstar(r_ind,:);qhstar(r_ind+Nk_in,:)];
    
    uhat = uhatstar(1:k_out+1);
    
end