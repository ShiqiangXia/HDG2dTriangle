function Plot2D(dom_type, x_pts,y_pts,f_vals,legend_str,save_flag, save_text)
    
    % Assume the local points are in increasing order
    
    [num_pts, num_ele] = size(x_pts);
    
    if strcmp(dom_type, 'Rec')
        N0 = sqrt(num_ele);
        % define data matrix
        total_pts = N0*num_pts;
        xpoint = zeros(total_pts,total_pts,numeric_t);
        ypoint = zeros(total_pts,total_pts,numeric_t);
        data_point = zeros(total_pts,total_pts,numeric_t);
        
        % go through each element
        % compute physical pts---> xpoint, ypoint
        for iy = 1:N0
            for ix = 1:N0
                ele_idx = (iy-1)*N0+ix;
%                 hx = element_list(ele_idx).len_x;
%                 hy = element_list(ele_idx).len_y;
%                 x_mid = element_list(ele_idx).x_mid;
%                 y_mid = element_list(ele_idx).y_mid;
                
                ix_star = (ix-1)*num_pts+1;
                ix_end = ix_star+num_pts-1;
                iy_star = (iy-1)*num_pts+1;
                iy_end = iy_star+num_pts-1;
                
                data_point(iy_star:iy_end,ix_star:ix_end) = f_vals(1:num_pts,1:num_pts,ele_idx);
                xpoint(iy_star:iy_end,ix_star:ix_end) =( (x_pts(:,ele_idx)).*ones(1,num_pts,numeric_t))';
                ypoint(iy_star:iy_end,ix_star:ix_end) = (y_pts(:,ele_idx)).*ones(1,num_pts,numeric_t);
            end
        end
        
        
        
        
    elseif strcmp(dom_type, 'Lshape')
        %{
        total_pts = 2*N0*num_pts;
        xpoint = zeros(total_pts,total_pts,numeric_t);
        ypoint = zeros(total_pts,total_pts,numeric_t);
        data_point = zeros(total_pts,total_pts,numeric_t)+nan;
        
       
        for iy = 1:2*N0
            for ix = 1:2*N0
                
                [flag,ele_idx]= CheckElementIndex(dom_type,iy,ix,N0);
               
                
                if flag == 1
                    hx = element_list(ele_idx).len_x;
                    hy = element_list(ele_idx).len_y;
                    x_mid = element_list(ele_idx).x_mid;
                    y_mid = element_list(ele_idx).y_mid;

                    ix_star = (ix-1)*num_pts+1;
                    ix_end = ix_star+num_pts-1;
                    iy_star = (iy-1)*num_pts+1;
                    iy_end = iy_star+num_pts-1;
                    
                    data_point(iy_star:iy_end,ix_star:ix_end) = f_vals(1:num_pts,1:num_pts,ele_idx)';
                    xpoint(iy_star:iy_end,ix_star:ix_end) =( (hx/numeric_t('2')*local_pts+x_mid).*ones(1,num_pts,numeric_t))';
                    ypoint(iy_star:iy_end,ix_star:ix_end) = (hy/numeric_t('2')*local_pts+y_mid).*ones(1,num_pts,numeric_t);
                    
                end  
            end
        end
        

     %}
        error('Lshape not implemented')
    else
        error("Mesh type not implemented")
        
    end
    
    xpoint = xpoint(end:-1:1,:);
    ypoint = ypoint(end:-1:1,:);
    data_point = data_point(end:-1:1,:);

     % surf plot
    p = figure;
    surf(xpoint,ypoint,data_point,'FaceAlpha',0.9);
    colormap jet
    shading interp
    xlabel("x");
    ylabel("y");
    title(legend_str,'Interpreter','latex','FontSize',18);
    colorbar;
    
    if save_flag == 1
        savefig(p,save_text);
    end

    
   
    
    
    
   
    
    
    
end