function [err_Jh,err_Jh_AC,err_Jh_elewise]= Error_Functional(func_type,para,...
                                    mymesh,GQ1DRef_pts,GQ1DRef_wts,...
                                    Jh,Jh_AC,Jh_elewise_list)
    NGQ = length(GQ1DRef_pts);
    num_elements = mymesh.num_elements;
    err_Jh_elewise = zeros(num_elements,1,numeric_t);
    
    if strcmp(func_type,'1')  % J(u) = (u,g)
        % need uexact, and source_g
        
        [uexact,source_g]=MyParaParse(para,'uexact','source_g');
        
        
        
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
            temp_int = Jk*GQ1DRef_wts'*(g_VD.*uexact_pts.*Jacobian_rs_to_ab)*GQ1DRef_wts;
            err_Jh_elewise(element_idx,1) = temp_int - Jh_elewise_list(element_idx,1);
            J_exact = J_exact + temp_int;
                  
            
        end
         % part 1: compute error J(u) - J(uh) and J(u) - Jh_AC 
         
         
              
    elseif strcmp(func_type,'2')  % J(u) = <-gradu * n,vD> = <q*n, vD>
        
        [qexact_1,qexact_2,vD]=MyParaParse(para,'qexact_1','qexact_2','vD');
        
        % go through all the faces, 
        % if it's a boundary face , do the integral
        
        J_exact = 0.0;
        
        for element_idx = 1: num_elements
            ele_face_idx_list  = mymesh.element_faces_list(element_idx,:);
            
            for ii = 1:length(ele_face_idx_list)
                face_id = ele_face_idx_list(ii);
                
                bdry_flag = mymesh.f_type(face_id);
                if bdry_flag == 1
                    temp_element = mymesh.element_list(element_idx,:);
                    V1 = mymesh.vertices_list(temp_element(ii),:);
                    if(ii==3)
                        V2_id = 1;
                    else
                        V2_id = ii+1;
                    end
                    V2 = mymesh.vertices_list(temp_element(V2_id),:);
                    
                    GQ_face_pts = 0.5*(V2-V1).*GQ1DRef_pts + 0.5*(V2+V1);
                    
                    q_face_pts1 = qexact_1(GQ_face_pts);
                    q_face_pts2 = qexact_2(GQ_face_pts);

                    vD_face_pts = vD(GQ_face_pts);

                    % normal verctor n?
                    Rot_mat = [0,-1;1,0]; % counterclock wise rotate 90 degree
                    normal_vec = - Rot_mat*(V2'-V1');
                    normal_vec = normal_vec/VectorNorm(normal_vec);

                    temp_formula = (q_face_pts1 * normal_vec(1)...
                        + q_face_pts2 * normal_vec(2)).*vD_face_pts;
                    
                    face_length = VectorNorm(V1'-V2');
                    
                    J_exact = J_exact  + 0.5 * face_length * GQ1DRef_wts' * temp_formula;
                 end
                
            end
            
        end
        
 
    else
        error ('This type of error has not been implemented yet.')
    end
    
    % J_exact = pi^2/2; % = int(sin(pi x)sin(pi y) * 2*pi^2* sin(pi x)sin(pi y))
    %J_exact = 4/(pi^2);  % = int(sin(pi x)sin(pi y) * 1)
    %J_exact = 4*(-4+pi^2)/(pi^4); % = int(sin(pi x)sin(pi y) * (x^2+y^2))
    
    %J_exact = 9/80 * (1+2*2^(1/3)); % = int (u_corner_singular * 1 )
    
    %J_exact = 3.173414780441447 ; % = int (u_corner_singular * 2*pi^2*sin(pi x)*sin(pi y) )
    %J_exact = 0.3088474886380344;  % = int (u_corner_singular * (x^2+y^2) )
    
    
    err_Jh = (J_exact - Jh);
    
    err_Jh_AC = (J_exact - Jh_AC);
    

end