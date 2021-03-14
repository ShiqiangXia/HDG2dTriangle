function [err_Jh,err_Jh_AC]= Error_Functional(func_type,para,...
                                    mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                                    Jh,Jh_AC,err_terms_flag)
    
    if strcmp(func_type,'1')  % J(u) = (u,g)
        % need uexact, and source_g
        
        [uexact,source_g]=MyParaParse(para,'uexact','source_g');
        
        NGQ = length(GQ1DRef_pts);
        num_elements = mymesh.num_elements;
        % Get Gauss Quadpoints on the square
        [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts);
        % Map Gauss Quadpoints to the reference triangle
        [r_list,s_list] = ABtoRS(a_list,b_list);
        
        J_exact = 0.0;
        % compute (u,g) over the mesh elementes
        for element_idx = 1: num_elements
            
            temp_element = mymesh.element_list(element_idx,:);
            vertice_list = mymesh.vertices_list(temp_element(:),:);
            Jk = mymesh.Jacobian_list(element_idx);
            % Gauss points on any element
            [x_list,y_list] = RStoXY(r_list,s_list,Jk,vertice_list);

            g_VD = source_g([x_list,y_list]);
            g_VD = reshape(g_VD,[],NGQ);
            
            uexact_pts = uexact([x_list,y_list]);
            uexact_pts = reshape(uexact_pts,[],NGQ);
            
            J_exact = J_exact +...
                Jk*GQ1DRef_wts'*(g_VD.*uexact_pts.*Jacobian_rs_to_ab)*GQ1DRef_wts;  
            
        end
         % part 1: compute error J(u) - J(uh) and J(u) - Jh_AC 
         
         err_Jh = abs(J_exact - Jh);
         err_Jh_AC = abs(J_exact - Jh_AC);
        
        if err_terms_flag == 1
            % part 2: compute error terms based on formula
            
            % need qexact, pexact
            [qexact_1,qexact_2,pexact_1,pexact_2]= MyParaParse(para,...
                'qexact_1','qexact_2','pexact_1','pexact_2');
            
            % compute the error terms over all elements
            
            % err_1 (q-qh,p-ph)
            % err_2 (q-qh,ph+grad(vh)) + (qh+grad(uh),p-ph)
            % err_3 <qhat*n - q*n, vh-vhat> + < uh- uhat, phat*n-p*n>

        end
        
        
    else
        error ('This type of error has not been implemented yet.')
    end
   
    
    
    
    
end