function or_list = GetUhatOritation(e)
    
    or_list = [e(:,2)>e(:,1), e(:,3)>e(:,2), e(:,1)>e(:,3)];
    
end