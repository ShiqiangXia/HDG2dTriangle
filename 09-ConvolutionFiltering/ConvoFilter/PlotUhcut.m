function PlotUhcut(uh_pts, hx,y_cut,N_level, GQ_x,legend_text,title_text )
    % only works for square so far
    
    [NGQ,num_ele] = size(GQ_x);
    xpts = reshape(GQ_x,[],1);
    xpts = xpts';
    fpts = zeros(1, num_ele*NGQ, numeric_t);
    
    for ii = 1:num_ele
        fpts(1,(ii-1)*NGQ+1:ii*NGQ) = uh_pts(y_cut,:,ii);
    end
    
    N_col = sqrt(num_ele);
    stard_id = (N_level-1)*N_col*NGQ + 1;
    end_id = N_level*N_col*NGQ;
    
    
    figure
    plot(xpts(stard_id:end_id),fpts(stard_id:end_id),'*--')
    hold on
    yline(0,'r-');
    for j = 1:N_col-1
        xline(hx*j,'r--');
    end
    
    title(title_text)
    legend(legend_text)
    
    
end