function res = get_order1_elements(col,row,Nx,Ny,N_corner_x,N_corner_y)
    % consider element (col(t), row(t)) as the center
    % expand based on N_corner_x and N_corner_y,  
    % namely: col(t)-N_corner_y, col(t) + N_corner_y,  of course we need to
    % consider the boundary Nx, Ny
    % record the covered elements
    m = length(col);
    res =[];
    for t = 1:m
        x_ind = row(t);
        y_ind = col(t);
        if x_ind - N_corner_x + 1 >= 1
            x_l = x_ind - N_corner_x + 1;
        else
            x_l = 1;
        end
        if x_ind + N_corner_x -1 <= Nx
            x_r = x_ind + N_corner_x -1;
        else
            x_r = Nx;
        end

        if y_ind - N_corner_y+1 >= 1
            y_l = y_ind - N_corner_y+1;
        else
            y_l = 1;
        end

        if y_ind + N_corner_y-1<=Ny
            y_r = y_ind + N_corner_y-1;
        else
            y_r = Ny;
        end
        
        nx = (x_r - x_l+1);
        ny = (y_r - y_l+1);
        Num = nx*ny;
        temp = zeros(Num,2,numeric_t);
        for j = y_l:y_r
            for i = x_l:x_r
                temp( (j-y_l)*nx+i-x_l+1,1) = j; % y index
                temp( (j-y_l)*nx+i-x_l+1,2) = i; % x index

            end
        end
        
        res = [res; temp];
    
        
    end
    
    res = unique(res,'rows');
    
       
    
end