function mymesh = build_mesh_by_mask(dom_type,coarse_mesh,mask,Nx_coarse,Ny_coarse)
    
    if strcmp(dom_type,'Rec')
        relation_mat = BuildRelationof2Meshes(coarse_mesh, coarse_mesh);
        element_list = [];
        Nsquare = Ny_coarse * Nx_coarse;
        for j = 1:Ny_coarse
            for i = 1:Nx_coarse
                if mask(j,i) == 1
                    square_idx = (j-1)*Nx_coarse + i;
                    temp_element_list = [relation_mat{square_idx},relation_mat{square_idx + Nsquare}];
                    for t = 1: length(temp_element_list)
                        ele_id = temp_element_list(t);
                        element_list = [element_list; coarse_mesh.element_list(ele_id,:)];
                    end
                end
            end
        end
        
        [p_new,e_new] = fixmesh(coarse_mesh.vertices_list,element_list );
        e_new = Countclockwise(p_new,e_new); % make sure a counterclockwise ordering of the vertices
        e_new = LongestEdgeFirst(p_new,e_new); % make sure the longest edge is the first

        bdry_edges = boundedges(p_new,e_new);
        % find boundary nodes
        tol = 1e-8;
        unit_square = @(p)drectangle(p,0,1,0,1);
        nodes_outer = find(abs(unit_square(p_new))<tol );% distance = 0  means on the boundary
        
        OUT= 1;
        INN = 11;
        
        [f,ef,f_type] = LabelFaces2(e_new, bdry_edges, nodes_outer,OUT,INN);
        new_dom_type = "Comp_Rec";
        mymesh = Mesh(new_dom_type, p_new, e_new, f, ef, f_type,"","");
        
        
    end

end