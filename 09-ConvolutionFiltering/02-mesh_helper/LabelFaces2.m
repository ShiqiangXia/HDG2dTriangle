function [f,ef,f_type] = LabelFaces2(e,...
        bdry_faces,nodes_bdry1,bdry_id1,bdry_id2)
    
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
%     bdry_id1 =;
%     bdry_id2 = ;
    
    for ii = 1:num_bd_f
        
        % the index of two vertices of one edge
        temp_f = bdry_faces(ii,:);
        
        
        bdry1_flag  = sum(ismember(temp_f,nodes_bdry1));
        %bdry2_flag  = sum(ismember(temp_f,nodes_bdry2));
        
        if bdry1_flag == 2 % this face is on the 1st boundary
            
            % find face index
            
            f_idx = find(ismember(f,temp_f,'rows'));
            f_type(f_idx) = bdry_id1;
            
%         elseif bdry2_flag == 2 % this face is on the 2nd boundary
        else
            f_idx = find(ismember(f,temp_f,'rows'));
            f_type(f_idx) = bdry_id2;
            
        end
        
    end
    
    

   
   
    
end