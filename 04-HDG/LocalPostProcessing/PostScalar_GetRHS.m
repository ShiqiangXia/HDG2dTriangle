function W_RHS = PostScalar_GetRHS(k,Jk,vertice_list,Auur,Auus,qh,uh)
    
    Nk = (k+1)*(k+2)/2;
    
    [~,Inv_AffineMap] = Ref_Tri_Map(Jk,vertice_list);
    r_x = Inv_AffineMap(1,1);
    r_y = Inv_AffineMap(1,2);
    s_x = Inv_AffineMap(2,1);
    s_y = Inv_AffineMap(2,2);
    
    r_ind = Hierarchical_indx(k+1);
    
    Buur = Auur(r_ind,:); 
    Buus = Auus(r_ind,:);
    
    Buru = Buur';
    Busu = Buus';
    
    W_RHS = Jk* (Buru*r_x+Busu*s_x) * qh(1:Nk)...
            + Jk * (Buru*r_y+Busu*s_y) * qh(Nk+1:end);
        
    W_RHS(1,1)=  uh(1)*Jk; 
    
    
end