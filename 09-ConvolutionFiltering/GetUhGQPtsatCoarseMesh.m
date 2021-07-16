function uh_pts = GetUhGQPtsatCoarseMesh(dom_type, k, coarse, fine, uh, GQ1DRef_pts)
    
    % Get relation mat
    relation_mat = BuildRelationof2Meshes(coarse, fine);
    
    num_coarse = coarse.num_elements ;
    Nsquare = num_coarse / 2; % two triangels make a square
    
    NGQ = length(GQ1DRef_pts);
    % store the uh values for each element
    uh_pts = zeros(NGQ, NGQ, Nsquare, numeric_t); 
    
    for ii = 1:Nsquare
        % go through each square
        
        % the triangels in the fine mesh 
        candidate_triangles = [relation_mat{ii}, relation_mat{ii + Nsquare}] ;
        
        % GQ points for this physical squaere
        [GQ_x, GQ_y] = GetPhyGQPts(dom_type, ii, Nsquare, GQ1DRef_pts);
        
        uh_pts(:,:,ii) = GetUhvales(k,uh,fine, GQ_x, GQ_y, candidate_triangles);
        
    end
    
    
end