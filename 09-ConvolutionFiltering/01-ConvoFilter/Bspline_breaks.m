function [breaks,break_flag] = Bspline_breaks(deg_bspine,center_pt)
% get spline breaks 
% see details in Convolution_mtraix.m

knots = (-deg_bspine:2:deg_bspine)*numeric_t('0.5');
s_l = center_pt - numeric_t('0.5');
s_r = center_pt + numeric_t('0.5');

if s_l >= knots(end) || s_r<=knots(1)
    break_flag = 1;
    breaks = 0;
    
elseif s_l< knots(1) && s_r>knots(1)
    break_flag = 2;
    breaks = [knots(1),s_r];
    
elseif s_r>knots(end) && s_l<knots(end)
    break_flag = 2;
    breaks = [s_l,knots(end)];
else
    temp_idx = find(knots>=s_l,1);
    temp_pt = knots(temp_idx);
    if (temp_pt == s_l)
        break_flag = 4;
        breaks = [s_l,s_r];
    else
        break_flag = 3;
        breaks = [s_l,temp_pt,s_r];
    end
end

end
