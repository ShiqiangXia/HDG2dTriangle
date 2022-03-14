function time_format(t_total, text)
    t_h = floor(t_total/3600);
    t_min = floor((t_total - t_h*3600)/60);
    t_sec = t_total - t_h*3600 - t_min*60 ;
    fprintf(text)
    fprintf(' %d h %d min %.1f s\n',t_h,t_min,t_sec);
    
end