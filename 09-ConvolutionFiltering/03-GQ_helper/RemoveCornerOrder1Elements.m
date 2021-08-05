function mat = RemoveCornerOrder1Elements(dom_type, mat,N_corner_x,N_corner_y,Nx, Ny)
    
    if  strcmp(dom_type,'Rec')
        
        for m1 = 1:1:N_corner_y % y-direction index
            for  m2 = 1:1:N_corner_x % x-direction index
                
                  mat(:,:,(m1-1)*Nx+m2) = 0;
                  
            end
        end
    end

end