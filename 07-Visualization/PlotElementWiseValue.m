function PlotElementWiseValue(mymesh,estimator,title_text)
    
    p = mymesh.vertices_list;
    e = mymesh.element_list;
    figure;
    myfigure = trisurf(e,p(:,1),p(:,2),0*p(:,1));
    temp = abs(estimator);
    max_err = max(temp);
    
    color_list = temp/max_err;
    
    set(gca,'CLim',[min(color_list), max(color_list)]);
    
    set(myfigure,'FaceColor','flat',...
       'FaceVertexCData',color_list,...
       'CDataMapping','scaled');
   colormap(jet);
   colorbar;
   view(2)
   axis equal
   ax=axis;axis(ax*1.001);
   title(title_text);

    
end