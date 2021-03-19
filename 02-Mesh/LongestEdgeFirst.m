function e = LongestEdgeFirst(p,e)
    % make sure the first edge is the longest edge of a triangle
    % (need this to make usre the refinement is longest edge based 
    % refinement and keep shape regularity)
    
    V1_list = p(e(:,1),:);
    V2_list = p(e(:,2),:);
    V3_list = p(e(:,3),:);
    
    edgelength(:,1) = sum((V1_list-V2_list).^2,2);
    edgelength(:,2) = sum((V2_list-V3_list).^2,2);
    edgelength(:,3) = sum((V3_list-V1_list).^2,2);
    
    [temp,I]=max(edgelength,[],2);
    
    e((I==2),[1 2 3])=e((I==2),[2 3 1]);
    e((I==3),[1 2 3])=e((I==3),[3 1 2]);
    
    
    
end