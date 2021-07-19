function out = Bspline_int(deg_bspine,center_pt,Nu,r,w)


% get breaks
[breaks,break_flag] = Bspline_breaks(deg_bspine,center_pt);

% integration for each break
if break_flag == 1
    out = numeric_t('0');
elseif break_flag == 2 ||break_flag==4
    x_l = breaks(1);
    x_r = breaks(2);
    % compute intetral
    out = Bspine_break_quadrature(deg_bspine,center_pt,x_l,x_r,r,w,Nu);
    
else
    x_l = breaks(1);
    x_m = breaks(2);
    x_r = breaks(3);
     % compute two intetrals
    out = Bspine_break_quadrature(deg_bspine,center_pt,x_l,x_m,r,w,Nu)+Bspine_break_quadrature(deg_bspine,center_pt,x_m,x_r,r,w,Nu);
    
end

end