function [f,ef,f_type] = LabelFaces(e,...
        nodes_dirichlet,nodes_neuman)
    
    % step 1: Label faces
    
    ne = size(e,1);
    % Get all edges
    enumerate = [1,2,2,3,3,1];
    faces = reshape(e(:,enumerate),[],2); % total 3*ne faces

    [f,~,ie] = unique(sort(faces,2),'rows'); % remove repeated faces

    % ie is the index of new labeled faces
    ef = reshape(ie(1:3*ne),[],3);
    
    % step 2: Mark boundaries.
     
    % go through all faces and record the boundary.
    
    num_f = length(f);
    f_type = zeros(num_f,1);
    
    INTERIOR = 0;
    DIRI = 1;
    NEU = 2;
    
    for ii = 1:num_f
        
        % the index of two vertices of one edge
        v1 = f(ii,1);
        v2 = f(ii,2);
        
        dirichlet_flag  = sum(ismember([v1 v2],nodes_dirichlet));
        neuman_flag  = sum(ismember([v1 v2],nodes_neuman));
        
        if dirichlet_flag == 2 % this face is on the boundary
            
            f_type(ii) = DIRI;
            
        elseif neuman_flag == 2
            
            f_type(ii) = NEU;
            
        else
            f_type(ii) = INTERIOR;
        end
        
    end
    
    

   
   
    
end