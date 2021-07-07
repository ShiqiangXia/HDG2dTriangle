function [err,err_list] = L2Error_qhs(mymesh, k1, qh1, k2, qh2, GQ1DRef_pts,GQ1DRef_wts, flag)
    % compute the error of qh1 - qh2 (they may have different degrees)
    
    if flag == 0
        % qh1 and qh2 are in Pk1 and Pk2
        Nu1 = (k1+2)*(k1+1)/2;
        Nu2 = (k2+2)*(k2+1)/2;

        qh1_1 = qh1(1:Nu1,:);
        qh1_2 = qh1(Nu1+1:end,:);

        qh2_1 = qh2(1:Nu2,:);
        qh2_2 = qh2(Nu2+1:end,:);

        [err1,err_list1] = L2Error_uhs(mymesh,k1, qh1_1, k2, qh2_1, ...
            GQ1DRef_pts,GQ1DRef_wts1,0) ;

        [err2,err_list2] = L2Error_uhs(mymesh,k1, qh1_2, k2, qh2_2,...
            GQ1DRef_pts,GQ1DRef_wts,0) ;

        
    elseif flag == 1
        % qh1 is in RTk1, qh2 is in Pk2
        Nu2 = (k2+2)*(k2+1)/2;
        qh2_1 = qh2(1:Nu2,:);
        qh2_2 = qh2(Nu2+1:end,:);
        
        [err1,err_list1] = L2Error_uhs(mymesh,k1, qh1, k2, qh2_1, ...
            GQ1DRef_pts,GQ1DRef_wts1,1);

        [err2,err_list2] = L2Error_uhs(mymesh,k1, qh1, k2, qh2_2,...
            GQ1DRef_pts,GQ1DRef_wts,2) ;
    else
        error('wrong flag')
    end
    
    err = sqrt((err1^2+err2^2));
    err_list = sqrt((err_list1.^2+err_list2.^2));
end