function uh_pts = GetGradUhPts(flag,Jk,vertice_list,V2D_r,V2D_s,coeff)
    
    % get grauh_x or graduh_y based on coeff of uh
    
    [~,Inv_AffineMap] = Ref_Tri_Map(Jk,vertice_list);
    r_x = Inv_AffineMap(1,1);
    r_y = Inv_AffineMap(1,2);
    s_x = Inv_AffineMap(2,1);
    s_y = Inv_AffineMap(2,2);
    
    if flag ==1  % u_x
        uh_pts = (V2D_r * r_x + V2D_s * s_x) * coeff;
    elseif flag == 2 % u_y
        uh_pts = (V2D_r * r_y + V2D_s * s_y) * coeff;
    else
        error('Wrong flag')
    end
    
    
    
end