function [err,err_list] = L2Error_scalar(mymesh,uh,...
                GQ1DRef_pts,GQ1DRef_wts,post_flag,...
                order,uexact)
            
            % compute the L2 error of ||uexact-uh||
            
            % get basic info
            NGQ = length(GQ1DRef_pts);
            num_elements = mymesh.num_elements;
            err_list = zeros(num_elements,1,numeric_t);
            
            k = order;
            [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts);

            % set up vandermonde matrix of uh
            if post_flag == 0
                V2D = Vandermonde2D(k,a_list,b_list);  % scalar Pk
            elseif post_flag == 1
                V2D = Vandermonde2D(k+1,a_list,b_list);  % scalar P_{k+1}
            elseif post_flag == 3
                [V2D,~] = RTVandermonde2D(k,a_list,b_list); % 1st comp. of RT_k
            elseif post_flag == 4
                [~,V2D] = RTVandermonde2D(k,a_list,b_list);% 2nd comp. of RT_k
            elseif post_flag == 5 || post_flag == 6
                [V2D_r, V2D_s] = GradVandermonde2D(k,a_list,b_list);
            else
                error('Wrong post-flag')
            end
            
            [r_list,s_list] = ABtoRS(a_list,b_list);
                      
            % go through each element, compute uexact at quadrature points
            for element_idx = 1: num_elements
                
                temp_element = mymesh.element_list(element_idx,:);
                vertice_list = mymesh.vertices_list(temp_element(:),:);
                Jk = mymesh.Jacobian_list(element_idx);
                
                [x_list,y_list] = RStoXY(r_list,s_list,Jk,vertice_list);
                
                uexact_pts = uexact([x_list,y_list]);
                
                uexact_pts = reshape(uexact_pts,[],NGQ);
                
                if post_flag < 5
                    uh_pts = V2D * (uh(:,element_idx));
                elseif post_flag == 5
                    uh_pts = GetGradUhPts(1,Jk,vertice_list,V2D_r,V2D_s,uh(:,element_idx));
                elseif post_flag == 6
                    uh_pts = GetGradUhPts(2,Jk,vertice_list,V2D_r,V2D_s,uh(:,element_idx));
                end
                
                uh_pts = reshape(uh_pts,[],NGQ);
                
                diff2 = (uexact_pts - uh_pts).^2;
                
                % do quadratrue.
                err_list(element_idx,1) =...
                    Jk*GQ1DRef_wts'*(diff2.*Jacobian_rs_to_ab)*GQ1DRef_wts;
  
            end
            
            err = sqrt(sum(err_list));
            err_list = sqrt(err_list);
            
           
end
