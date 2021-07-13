function idx = FindCoarseElement(pts,num, element_list, vertices_list)
    
    idx = 0;
    
    for ii = 1:num
        P1 = vertices_list(element_list(ii,1),:);
        P2 = vertices_list(element_list(ii,2),:);
        P3 = vertices_list(element_list(ii,3),:);
        
        flag_3pts = true;
        
        for jj = 1:3
            % check if all three points are inside the triangle
            flag = CheckPointInTriangle(pts(jj,:), P1,P2,P3);
            if flag
                continue;
            else
                flag_3pts = false;
                break;
            end
        end
        
        if flag_3pts
            idx = ii;
            break
        end
    end
end