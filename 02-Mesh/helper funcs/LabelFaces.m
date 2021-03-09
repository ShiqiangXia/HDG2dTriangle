function [f,ef,f_type] = LabelFaces(e,...
        nodes_dirichlet,nodes_neuman,bdry_faces)
    
    % step 1: Label faces
    
    ne = size(e,1);
    % Get all edges
    enumerate = [1,2,2,3,3,1];
    faces = reshape(e(:,enumerate)',2,[]); % total 3*ne faces
    faces = faces';

    [f,~,ie] = unique(sort(faces,2),'rows'); % remove repeated faces

    % ie is the index of new labeled faces
    ef = reshape(ie(1:3*ne),3,[]);
    ef=ef';
    
    % step 2: Mark boundaries.
     
    % go through all faces and record the boundary.
    
    bdry_faces = sort(bdry_faces,2);
    num_bd_f = length(bdry_faces);
    
    num_f = length(f);
    f_type = zeros(num_f,1);
    
    %INTERIOR = 0;
    DIRI = 1;
    NEU = 2;
    
    for ii = 1:num_bd_f
        
        % the index of two vertices of one edge
        temp_f = bdry_faces(ii,:);
        
        
        dirichlet_flag  = sum(ismember(temp_f,nodes_dirichlet));
        neuman_flag  = sum(ismember(temp_f,nodes_neuman));
        
        if dirichlet_flag == 2 % this face is on dirichlet the boundary
            
            % find face index
            
            f_idx = find(ismember(f,temp_f,'rows'));
            f_type(f_idx) = DIRI;
            
        elseif neuman_flag == 2 % this face is on neuman the boundary
            
            f_idx = find(ismember(f,temp_f,'rows'));
            f_type(f_idx) = NEU;
            
        end
        
    end
    
    

   
   
    
end