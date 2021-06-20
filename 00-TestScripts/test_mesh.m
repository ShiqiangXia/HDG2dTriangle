global class_t;
class_t = "double";

% rectangle mesh

mymesh = Build2DMesh(1,'Rec',0.5,["bottom","top","left","right"],[],0,0,5,5,0);

mymesh.Plot(1);

% uniform refine
mesh2 = mymesh.UniformRefine();
mesh2.Plot(1);

% only refine selected elements
mesh3 = mymesh.Refine([1,2,3],'RGB');
mesh3.Plot(0);

% L shpaed mesh
mymesh = Build2DMesh(1,'L',0.2,["bottom","top_high","right_low","left"],["right_high","top_low"]);
mymesh.Plot(1)

% Cicle mesh

%mymesh = Build2DMesh(0,'Cir',0.3,["Circle"],[],0,0,1);

% Cicle with hole
%mymesh = Build2DMesh(0,'CirHole',0.1,["InnerCircle"],["OuterCircle"],0,0,1,0,0,0.5);