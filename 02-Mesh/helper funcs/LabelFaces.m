function [f,ef,f_type] = LabelFaces(dom_type,p,e,...
        dirichlet_flag,neuman_flag,varargin)
    
    % step 1: Label faces
    
    ne = size(e,1);
    % Get all edges
    enumerate = [1,2,2,3,3,1];
    faces = reshape(e(:,enumerate),[],2); % total 3*ne faces

    [f,~,ie] = unique(sort(faces,2),'rows'); % remove repeated faces

    % ie is the index of new labeled faces
    ef = reshape(ie(1:3*ne),[],3);
    
    % step 2: Mark boundaries.
    
    tol = 1e-8;
    
    % find boundary nodes
    if isempty(intersect(dirichlet_flag,neuman_flag))~=1
        error("The dirichlet and neuman boudary should not intersect!")
        
    end
    
    fd_dirichlet = GetDistFunctions(dom_type,dirichlet_flag,varargin{:});
    nodes_dirichlet = find(abs(fd_dirichlet(p))<tol );% distance = 0  means on the boundary
    
    fd_neuman = GetDistFunctions(dom_type,neuman_flag,varargin{:});
    nodes_neuman = find(abs(fd_neuman(p))<tol );
    
    figure
    labs = 1:length(p);
    plot(p(:,1),p(:,2),'k.',...
        p(nodes_dirichlet,1),p(nodes_dirichlet,2),'rx',...
        p(nodes_neuman,1),p(nodes_neuman,2),'bd')
    legend('vertices','dirichlet', 'neuman')
    labelpoints(p(:,1),p(:,2),labs,'NW',0.5,0,'FontSize', 14);
    
    
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