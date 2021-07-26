function uh_pts = GetUhGQPtsatCoarseMesh(k, coarse, fine, uh, GQ_x, GQ_y)
    % GOAL: evaluate uh at GQ points defined on coarse mesh (square element)
    % here uh is defined on fine mesh (triangle element)
    % OUTPUT: Matrix of NGQ x NGQ x num_element
    % each row is the same y point different x points
 
    % all the GQ physical points are stored in GQ_x and GQ_y (NGQ x num_ele)
   
    % Get relation mat
    relation_mat = BuildRelationof2Meshes(coarse, fine);
    
    num_coarse = coarse.num_elements ;
    Nsquare = num_coarse / 2; % two triangels make a square
    
    [NGQ,~] = size(GQ_x);
    % store the uh values for each element
    uh_pts = zeros(NGQ, NGQ, Nsquare, numeric_t); 
    
    for ii = 1:Nsquare
        % go through each square
        
        % the triangels in the fine mesh 
        candidate_triangles = [relation_mat{ii}, relation_mat{ii + Nsquare}] ;
        
        uh_pts(:,:,ii) = GetUhvales(k,uh,fine, GQ_x(:,ii), GQ_y(:,ii), candidate_triangles);
        
    end
    
    
end