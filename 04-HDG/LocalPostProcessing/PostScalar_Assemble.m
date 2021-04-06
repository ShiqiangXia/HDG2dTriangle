function W_mat = PostScalar_Assemble(Jk,vertice_list,Aurur,Aurus,Ausus)
    
    
    [~,Inv_AffineMap] = Ref_Tri_Map(Jk,vertice_list);
    r_x = Inv_AffineMap(1,1);
    r_y = Inv_AffineMap(1,2);
    s_x = Inv_AffineMap(2,1);
    s_y = Inv_AffineMap(2,2);
   
    
    W_mat = Jk*(Aurur*r_x*r_x + Aurus*r_x*s_x + Aurus'*s_x*r_x +Ausus*s_x*s_x...
                   +Aurur*r_y*r_y + Aurus*r_y*s_y +Aurus'*s_y*r_y + Ausus*s_y*s_y  );
         
    
    W_mat(1,1) = Jk;
    
end