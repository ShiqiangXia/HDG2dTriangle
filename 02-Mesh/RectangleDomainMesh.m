function mymesh = RectangleDomainMesh(structured_flag,...
        x1,y1,x2,y2,h0,...
        dirichlet_flag,neuman_flag)
    
    % This code generates a triangle mesh for a rectange with corners
    % (x1,y1),(x2,y2) with size h
    % structured_flag: 0/1 structured mesh or not 
    % specify the bounary data type by 'bottom','left','right','top'
    
    if structured_flag == 0
        % unstructured mesh
        
        % step 1: Use distmesh to generatre unstructured meshes
        
        fd = @(p)drectangle(p,x1,x2,y1,y2);
        bbox = [x1,y1; x2,y2];
        pfix = [x1,y1;x2,y1;x2,y2;x1,y2];
        [p,e] = distmesh2d(fd, @huniform,h0,bbox,pfix);
        
        e = Countclockwise(p,e); % make sure a counterclockwise ordering of the vertices
        
        [f,ef,dirichlet,neuman,BCType] = RecLabelFaces(p,e,x1,x2,y1,y2,dirichlet_flag,neuman_flag);
        
        mymesh = Mesh('Rec',p,e,f,ef,dirichlet,neuman,BCType);
        
        
    elseif structured_flag == 1
        % structured mesh
        
    else
        error('Wrong structured_flag!')
        
    end
    
end
