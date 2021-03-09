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
            
            if nargin-5 < 4
                error("Not enought parameter to generate a mesh for Rectangle domain")
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
            
        elseif strcmp(dom_type,'Cir')
            
            if nargin -5 < 3
                error("Not enought parameter to generate a mesh for a circle")
            end
            
            x1=varargin{1};
            y1=varargin{2};
            r =varargin{3};
            
            fd = @(p)dcircle(p,x1,y1,r);
            
            bbox = [x1-r,y1-r; x1+r,y1+r];
            [p,e]=distmesh2d(fd,@huniform,h0,bbox,[]);
            
        elseif strcmp(dom_type,'CirHole')
            
            if nargin - 5 < 6
                error("Not enought parameter to generate a mesh for a circle")
            end
            
            x1=varargin{1};
            y1=varargin{2};
            r1 =varargin{3};
            
            x2=varargin{4};
            y2=varargin{5};
            r2 =varargin{6};
            
            if (x1+r1<=x2+r2 || y1+r1<= y2+r2)
                error("The inner cicle should be inside the outer circle. Check papamters! ")
            end
            
            fd = @(p)ddiff( dcircle(p,x1,y1,r1),dcircle(p,x2,y2,r2) );
            
            bbox = [x1-r1,y1-r1; x1+r1,y1+r1];
            
            [p,e]=distmesh2d(fd,@huniform,h0,bbox,[]);
            
            
        else
           
           error("Mesh of this domain type is not implemented yet!");
            
            
        end
        
    elseif structured_flag == 1
        % structured mesh
        
           if strcmp(dom_type,'Rec')
               if nargin-5 < 4
                    error("Not enought parameter to generate a mesh for Rectangle domain")
               end
            
               [p,e] = Structured2DMesh(dom_type,h0,varargin{:});
               simpplot(p,e);
               
           elseif strcmp(dom_type,'L')
               
               [p,e] = Structured2DMesh(dom_type,h0,varargin{:});
               simpplot(p,e);
               
           else
               error("This type of structured mesh is not implemented yet")
           end
        
    else
        error('Wrong structured_flag!')
        
    end
    
    
    e = Countclockwise(p,e); % make sure a counterclockwise ordering of the vertices
    
    % make sure dirichlet boundary and neuman boundary are different
    if isempty(dirichlet_flag) || isempty(neuman_flag)
        % do nothing
    elseif isempty(intersect(dirichlet_flag,neuman_flag))~=1
        error("The dirichlet and neuman boudary should not intersect!")
    end
    
    % find boundary nodes
    tol = 1e-8;
    if ~isempty(dirichlet_flag)
        fd_dirichlet = GetDistFunctions(dom_type,dirichlet_flag,varargin{:});
        nodes_dirichlet = find(abs(fd_dirichlet(p))<tol );% distance = 0  means on the boundary
    else
        nodes_dirichlet=[];
    end
    
    if ~isempty(neuman_flag)
        fd_neuman = GetDistFunctions(dom_type,neuman_flag,varargin{:});
        nodes_neuman = find(abs(fd_neuman(p))<tol );
    else
        nodes_neuman = [];
    end
    
    bdry_edges = boundedges(p,e);
    [f,ef,f_type] = LabelFaces(e,nodes_dirichlet,nodes_neuman,bdry_edges);

    mymesh = Mesh(dom_type,p,e,f,ef,f_type,dirichlet_flag,neuman_flag,varargin{:});
    
end
