function [err,err_list] = L2Error_vector(mymesh,qh,...
                GQ1DRef_pts,GQ1DRef_wts,post_flag,...
                order,qexact_1,qexact_2)
            
            k = order;
            Nu = (k+2)*(k+1)/2;
            
            qh_1 = qh(1:Nu,:);
            qh_2 = qh(Nu+1:end,:);
            
            % call scalar L2 error 
            if post_flag == 0
                [err1,err_list1] = L2Error_scalar(mymesh,qh_1,...
                    GQ1DRef_pts,GQ1DRef_wts,post_flag,...
                    order,qexact_1) ;

                [err2,err_list2] = L2Error_scalar(mymesh,qh_2,...
                    GQ1DRef_pts,GQ1DRef_wts,post_flag,...
                    order,qexact_2) ;

                err = sqrt((err1^2+err2^2));
                err_list = sqrt((err_list1.^2+err_list2.^2));
            elseif post_flag == 1
            else
                error('Error calculation for this type of post-processing has not implemented. ')
            end
            
            
end

            