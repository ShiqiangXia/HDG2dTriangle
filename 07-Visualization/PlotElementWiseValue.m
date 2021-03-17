function PlotElementWiseValue(mymesh,varargin)
    
    nn = (nargin - 1)/2;
    p = mymesh.vertices_list;
    e = mymesh.element_list;
    color_list = abs(varargin{1});
    
    for ii = 1:nn
        
        %subplot(nn,1,ii);
        figure;
        myfigure = trisurf(e,p(:,1),p(:,2),0*p(:,1));
        set(gca,'CLim',[min(color_list), max(color_list)]);
        temp  = varargin{2*ii-1};
        set(myfigure,'FaceColor','flat',...
           'FaceVertexCData',temp,...
           'CDataMapping','scaled');
       colormap(jet);
       colorbar;
       view(2)
       axis equal
       ax=axis;axis(ax*1.001);
       title(varargin{2*ii});
       
        
        
    end
    

%     myfigure = trisurf(e,p(:,1),p(:,2),0*p(:,1));
%     temp = abs(values_ele);
%     %temp = estimator;
%     %max_err = max(temp);
%     
%     color_list = temp;%/max_err;
%     
%     set(gca,'CLim',[min(color_list), max(color_list)]);
%     
%     set(myfigure,'FaceColor','flat',...
%        'FaceVertexCData',color_list,...
%        'CDataMapping','scaled');
%    colormap(jet);
%    colorbar;
%    view(2)
%    axis equal
%    ax=axis;axis(ax*1.001);
%    title(title_text);

    
end