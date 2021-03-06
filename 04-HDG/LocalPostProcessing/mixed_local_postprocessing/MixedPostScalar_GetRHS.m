function W_RHS = MixedPostScalar_GetRHS(Jk,vertice_list,RT1_ur, RT1_us, RT2_ur, RT2_us,qh_RT_coef, uh_ave)
    % build right hand side vector
    % each row i: sum : qh_j * (RT1_j, RT_j) * (ux_i, uy_i)
    % ux = ur * rx + us * sx
    % uy = ur * ry + us * sy
    
    [~,Inv_AffineMap] = Ref_Tri_Map(Jk,vertice_list);
    r_x = Inv_AffineMap(1,1);
    r_y = Inv_AffineMap(1,2);
    s_x = Inv_AffineMap(2,1);
    s_y = Inv_AffineMap(2,2);
    
    
    % RT1 * ux
    W_RHS_1 = Jk * (RT1_ur' * r_x + RT1_us' * s_x) * (-qh_RT_coef);
    % RT2 * uy
    W_RHS_2 = Jk * (RT2_ur' * r_y + RT2_us' * s_y) * (-qh_RT_coef);
    
    W_RHS = W_RHS_1 + W_RHS_2;
    W_RHS(1) = uh_ave*Jk; 
    
    
end