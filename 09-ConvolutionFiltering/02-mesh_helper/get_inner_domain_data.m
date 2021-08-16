function mat = get_inner_domain_data(dom_type, mat,Nx, Ny, mask)
    
    if  strcmp(dom_type,'Rec')
        for m1 = 1:1:Ny % y-direction index
            for  m2 = 1:1:Nx % x-direction index
                if mask(m1,m2) == 1
                    mat(:,:,(m1-1)*Nx+m2) = 0;
                end
            end
        end
    end
    
end