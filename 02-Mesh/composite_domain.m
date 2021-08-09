
x1 = 0;
y1 = 0;
x2 = 1;
y2 = 1;

unit_square = @(p)drectangle(p,x1,x2,y1,y2);


h = 0.05;
k = 2;
d = 0.1;
N_bd = 4;

x1_inner = x1+N_bd*h;
y1_inner = y1+N_bd*h;
x2_inner = x2-N_bd*h;
y2_inner = y2-N_bd*h;

inner_square = @(p)drectangle(p,x1_inner,x2_inner,y1_inner,y2_inner);

x1_corner = x1+N_bd*h;
y1_corner = y1+N_bd*h;
x2_corner = x1+d;
y2_corner = y1+d;

if x2_corner > x1_corner && y2_corner > y1_corner
    
    corner_square = @(p)drectangle(p,x1_corner,x2_corner,y1_corner,y2_corner);
    inner_domain = @(p)ddiff(inner_square(p),corner_square(p));
    
    bbox_inner = [x1_inner,y1_inner; x2_inner,y2_inner];
    pfix_inner = [x2_corner,y2_corner; x1_inner,y2_corner; ...
                 x1_inner,y2_inner; x2_inner,y2_inner;x2_inner,y1_inner;...
                 x2_corner,y1_inner];
    %[p_inner,e_inner] = distmesh2d(inner_domain, @huniform,h,bbox_inner,pfix_inner);
    
else
    inner_domain = inner_square;
    bbox_inner = [x1_inner,y1_inner; x2_inner,y2_inner];
    pfix_inner = [x1_inner,y1_inner;x2_inner,y1_inner;x2_inner,y2_inner;x1_inner,y2_inner];
end


bbox_outer = [x1,y1; x2,y2];
pfix_outer = [x1,y1;x2,y1;x2,y2;x1,y2;pfix_inner];
outer_domain = @(p)ddiff(unit_square(p),inner_domain(p));
[p_outer,e_outer] = distmesh2d(outer_domain, @huniform,h,bbox_outer,pfix_outer);

e_outer = Countclockwise(p_outer,e_outer); % make sure a counterclockwise ordering of the vertices
e_outer = LongestEdgeFirst(p_outer,e_outer); % make sure the longest edge is the first

bdry_edges = boundedges(p_outer,e_outer);

% find boundary nodes
tol = 1e-8;
nodes_outer = find(abs(unit_square(p))<tol );% distance = 0  means on the boundary
nodes_inner = find(abs(inner_domain(p))<tol );



