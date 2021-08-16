function mask = build_mesh_mask(dom_type, Nx, Ny, N_bd, special_idx)
    if strcmp(dom_type, 'Rec')
        mask = ones(Nx,Ny,numeric_t);
        
        for i = 1:Ny
            for j = 1:Nx
                if i>=N_bd+1 && i <= Ny-N_bd && j>=N_bd+1 && j<= Nx-N_bd
                    mask(i,j) = 0;
                end
            end
        end
        [n,~] = size(special_idx);
        for k=1:n
            mask(special_idx(k,1), special_idx(k,2)) = 1;
        end
    end
end