function relation_mat = BuildRelationof2Meshes(coarse, fine)
    
    % for each element in the fine mesh,
    % find what corase mesh it is contained in.
    % if the coarse mesh elements can't contain all fine mesh elements
    % return error
    
    num_coarse = coarse.num_elements;
    num_fine = fine.num_elements;
    relation_mat = cell(num_coarse,1);
    
    for ii = 1:num_fine
        temp_element = fine.element_list(ii,:);
        pts = fine.vertices_list(temp_element(:),:);
        
        % find which corase element is this belong to 
        
        idx = FindCoarseElement(pts, num_coarse, coarse.element_list, coarse.vertices_list);
        if idx == 0
            error('the two meshes are not subset relation')
        else
            relation_mat{idx,1} = [relation_mat{idx,1},ii];
        end
        
    end
    
end