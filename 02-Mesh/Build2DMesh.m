function mymesh = Build2DMesh(structured_flag,dom_type,...
        h0,...
        dirichlet_flag,neuman_flag,varargin)
   
    % This code generates a triangle mesh for different dom_type with size h
    %
    % dom_type = 'Rec', rectangle domain with (x1,y1),(x2,y2) 
    % dom_type = 'L', Lshaped domain[0,2]^2\[1,1]x[2,2], 
    % dom_type = 'Cir', circle at x1,y1 with radius r
    % 
    % varargin provides the parameters needed to define the domain
    %
    % dirichlet_flag, neuman_flag defined the data type of the boundary
    %
    % structured_flag: 0/1 structured mesh or not 
    % specify the bounary data type by 'bottom','left','right','top'
    
    if structured_flag == 0
        % unstructured mesh
        
        % step 1: Use distmesh to generatre unstructured meshes
        if strcmp(dom_type,'Rec')
            
            if nargin < 4
                error("Not enought parameter to generate mesh for Rectangle domain")
            end
            
            x1=varargin{1};
            y1=varargin{2};
            x2=varargin{3};
            y2=varargin{4};
            
            fd = @(p)drectangle(p,x1,x2,y1,y2);
            bbox = [x1,y1; x2,y2];
            pfix = [x1,y1;x2,y1;x2,y2;x1,y2];
            [p,e] = distmesh2d(fd, @huniform,h0,bbox,pfix);
            
        elseif strcmp(dom_type,'L')
            
            fd = @(p)ddiff(drectangle(p,0,2,0,2),drectangle(p,1,2,1,2));
            bbox = [0,0; 2,2];
            pfix = [0,0; 2,0; 2,1; 1,1; 1,2; 0,2];
            [p,e] = distmesh2d(fd, @huniform,h0,bbox,pfix);
            
        else
            
        end
        
        
        e = Countclockwise(p,e); % make sure a counterclockwise ordering of the vertices
        
        [f,ef,f_type] = LabelFaces(dom_type,p,e,dirichlet_flag,neuman_flag,varargin{:});
        
        mymesh = Mesh(dom_type,p,e,f,ef,f_type);
        
        
    elseif structured_flag == 1
        % structured mesh
        
    else
        error('Wrong structured_flag!')
        
    end
    
end
