function [x,y] = RStoXY(r_list,s_list,Jk,vertice_list)
    
    % [x,y]' = a_mat *[r+1,s+1]'+[x1,y1]'
    
    [AffineMap,~] = Ref_Tri_Map(Jk,vertice_list);
    
    V1 = vertice_list(1,:)';
    
    r_v = [r_list';s_list'] + [1;1];
    
    temp = AffineMap*r_v + V1;
    
    x = temp(1,:)';
    y = temp(2,:)';
    
end