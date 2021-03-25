function marked = ACh_ErrEstimate(ACh_elewise_list,tol_adp,percent, mark_flag)
    if mark_flag == 0
        max_err = max(abs(ACh_elewise_list));

        marked = abs(ACh_elewise_list) > percent*max_err ;
    elseif mark_flag == 1 % Dorfler marking strategy
        
        temp_err_est = abs(ACh_elewise_list);
        total_err_est = sum(temp_err_est);
        criterion = percent* total_err_est;
        
        [s,i] = sort(temp_err_est,'descend');
        nn = size(temp_err_est,1);
        marked = zeros(nn,1);
        part_sum = 0;
        
        for jj = 1:nn
            part_sum = part_sum + s(jj);
            if part_sum < criterion
                marked(i(jj)) = 1;
            else
                break;
            end
        end
        marked = logical(marked);
    else
        error('Marking strategy %d is not implemented.', mark_flag)
    end
    
end