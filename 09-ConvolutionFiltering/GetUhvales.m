function mat = GetUhvales(k,uh,mesh, GQ_x, GQ_y, candidate_triangles)
    
    n = lenghth(candidate_triangles);
    [NGQ,~] = size(GQ_x);
    mat = zeros(NGQ,NGQ, numeric_t);
    for jj = 1:NGQ
        for ii = 1:NGQ
            
            P = [GQ_x(ii),GQ_y(ii)];
            
            % check which triangle does point P belong to
            for s = 1:n
                ele_idx = candidate_triangles(s);
                P1 = mesh.vertices_list(mesh.element_list(ele_idx,1),:);
                P2 = mesh.vertices_list(mesh.element_list(ele_idx,2),:);
                P3 = mesh.vertices_list(mesh.element_list(ele_idx,3),:);
                
                flag = CheckPointInTriangle(P, P1,P2,P3);
                
                if flag
                    target_idx = s;
                    break;
                end
                
            end
            
            % evaluate uh in element target_idx at point P
            vertice_list = mesh.vertices_list(mesh.element_list(target_idx,:),:);
            Jk = mesh.Jacobian_list(target_idx);
            
            [r_list,s_list] = XYtoRS(P(1),P(2),Jk,vertice_list);
            [a_list,b_list] = RStoAB(r_list,s_list);
            V2D = Vandermonde2D(k,a_list,b_list);
            mat(jj,ii) = V2D * (uh(:,target_idx));
            
        end
    end
    
end