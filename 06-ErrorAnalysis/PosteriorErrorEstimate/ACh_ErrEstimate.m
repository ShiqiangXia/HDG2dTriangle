function marked = ACh_ErrEstimate(ACh_elewise_list,tol_adp,percent)
    
    max_err = max(abs(ACh_elewise_list));
    
    marked = abs(ACh_elewise_list) > percent*max_err ;
    
end