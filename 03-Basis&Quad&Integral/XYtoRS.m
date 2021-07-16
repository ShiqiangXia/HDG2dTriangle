function [r,s] = XYtoRS(x_list,y_list,Jk,vertice_list)
    
    % map point(x,y) from phy element to point (r,s ) on ref element
    % [r,s]' = a_mat^-1 *[x-x1,y-y1]'+[-1,-1]'
    
    
    [~,Inv_AffineMap] = Ref_Tri_Map(Jk,vertice_list);
    V1 = vertice_list(1,:)';
    
    pts_v = [x_list';y_list'] - V1;
    
    temp = Inv_AffineMap*pts_v + [-1;-1];
    
    r = temp(1,:)';
    s = temp(2,:)';
    
end