function [err,err_list] = L2Error_vector(mymesh,qh_coeff,...
                GQ1DRef_pts,GQ1DRef_wts,post_flag,...
                order,qexact_1,qexact_2)
            
            k = order;
            
            
            % call scalar L2 error 
            if post_flag == 0
                % since qh is (Pk,0) and (0,Pk)
                % we can split it as two scalar L2 errors
                % namely calculating qexact_1 - qh_1 and qexact_2 - qh_2
                Nu = (k+2)*(k+1)/2;
            
                qh_1 = qh_coeff(1:Nu,:);
                qh_2 = qh_coeff(Nu+1:end,:);
                
                [err1,err_list1] = L2Error_scalar(mymesh,qh_1,...
                    GQ1DRef_pts,GQ1DRef_wts,post_flag,...
                    order,qexact_1) ;

                [err2,err_list2] = L2Error_scalar(mymesh,qh_2,...
                    GQ1DRef_pts,GQ1DRef_wts,post_flag,...
                    order,qexact_2) ;

                
            elseif post_flag == 1
                % for this case. qh is in the Raviart-Thomas space Pk_vec + x*Pk_homo
                % unlike the previous case, this time qh_1 and qh_2 depends
                % on all the coefficents of qh and to evalue qh at
                % quadrature points correctly, we need to use the
                % Raviart-Thomas vandermonde matrix. This is achieved by
                % different calling L2Error_scalar with different post_flag
                
                [err1,err_list1] = L2Error_scalar(mymesh,qh_coeff,...
                    GQ1DRef_pts,GQ1DRef_wts,post_flag+2,...
                    order,qexact_1) ;
                
                [err2,err_list2] = L2Error_scalar(mymesh,qh_coeff,...
                    GQ1DRef_pts,GQ1DRef_wts,post_flag+3,...
                    order,qexact_2) ;
                
            elseif post_flag == 2
                % in this case the vector is -grad_uh = (-uh_x, -uh_y)
                
                [err1,err_list1] = L2Error_scalar(mymesh,qh_coeff,...
                    GQ1DRef_pts,GQ1DRef_wts,5,...
                    order,qexact_1) ;
                
                [err2,err_list2] = L2Error_scalar(mymesh,qh_coeff,...
                    GQ1DRef_pts,GQ1DRef_wts,6,...
                    order,qexact_2) ;
                
            else
                error('Error calculation for this type of post-processing has not implemented. ')
            end
            
            err = sqrt((err1^2+err2^2));
            err_list = sqrt((err_list1.^2+err_list2.^2));
            
            
end

            