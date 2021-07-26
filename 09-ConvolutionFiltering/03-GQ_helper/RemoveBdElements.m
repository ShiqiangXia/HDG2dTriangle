function mat = RemoveBdElements(dom_type, mat,N_bd,Nx, Ny)
    
    if  strcmp(dom_type,'Rec')
        
        for m1 = 1:1:Ny % y-direction index
            for  m2 = 1:1:Nx % x-direction index
                if m1 >= N_bd+1 && m1 <= Ny-N_bd && m2 >= N_bd+1 && m2 <= Nx-N_bd
                    continue
                else
                    mat(:,:,(m1-1)*Nx+m2) = 0;
                end
            end
        end
    end

end