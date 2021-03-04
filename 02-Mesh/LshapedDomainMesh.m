function mymesh = LshapedDomainMesh(structured_flag,...
        h0,...
        dirichlet_flag,neuman_flag)
    
    % This code generates a triangle mesh for a Lshaped domain [0,2]^2\[1,1]x[2,2]
    % with size h
    % structured_flag: 0/1 structured mesh or not 
    % specify the bounary data type by 'bottom','left','right','top'
    
    if structured_flag == 0
        % unstructured mesh
        
        % step 1: Use distmesh to generatre unstructured meshes
        
        fd = @(p)ddiff(drectangle(p,0,2,0,2),drectangle(p,1,2,1,2));
        bbox = [0,0; 2,2];
        pfix = [0,0; 2,0; 2,1; 1,1; 1,2; 0,2];
        [p,e] = distmesh2d(fd, @huniform,h0,bbox,pfix);
        
        e = Countclockwise(p,e); % make sure a counterclockwise ordering of the vertices
        
        [f,ef,f_type] = LshapeLabelFaces(p,e,dirichlet_flag,neuman_flag);
        
        mymesh = Mesh('L',p,e,f,ef,f_type);
        
        
    elseif structured_flag == 1
        % structured mesh
        
    else
        error('Wrong structured_flag!')
        
    end
    
end
