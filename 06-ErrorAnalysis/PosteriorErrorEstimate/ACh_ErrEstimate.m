function marked = ACh_ErrEstimate(Est_elewise_list,TOL,percent, mark_flag)
    if mark_flag == 0
        max_err = max(abs(Est_elewise_list));

        marked = abs(Est_elewise_list) > percent*max_err ;
        
    elseif mark_flag == 1 % Dorfler marking strategy
        
        temp_err_est = abs(Est_elewise_list);
        total_err_est = sum(temp_err_est);
        criterion = percent* total_err_est;
        
        [s,i] = sort(temp_err_est,'descend');
        nn = size(temp_err_est,1);
        marked = zeros(nn,1);
        part_sum = 0;
        
        for jj = 1:nn
            
            if part_sum < criterion
                marked(i(jj)) = 1;
            else
                break;
            end
            part_sum = part_sum + s(jj);
        end
        
        PercentPlot(s,jj);
        
        
    elseif mark_flag == 2 % equi distribute error
        
        temp_err_est = abs(Est_elewise_list);
        nn = size(temp_err_est,1);
        ele_tol = TOL/nn;
        marked = zeros(nn,1);
        
        for jj = 1:nn
            if temp_err_est(jj) > ele_tol
                marked(jj) = 1;
            end
        end
        
    elseif mark_flag == 3 % fixed fraction strategy
        
        temp_err_est = abs(Est_elewise_list);
        [s,i] = sort(temp_err_est,'descend');
        nn = size(temp_err_est,1);
        marked = zeros(nn,1);
        N_mark = ceil(percent*nn);
        for jj = 1:N_mark
            marked(i(jj)) = 1;
        end
        
        
    else
        error('Marking strategy %d is not implemented.', mark_flag)
    end
    
    marked = logical(marked);
    
end