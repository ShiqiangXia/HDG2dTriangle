
% rectangle mesh
mymesh = Build2DMesh(1,'Rec',1,["bottom","top","left","right"],[],0,0,1,1,0);

% L shpaed mesh
%mymesh = Build2DMesh(1,'L',0.2,["bottom","top_high","right_low","left"],["right_high","top_low"]);


% Cicle mesh

%mymesh = Build2DMesh(0,'Cir',0.3,["Circle"],[],0,0,1);

% Cicle with hole
%mymesh = Build2DMesh(0,'CirHole',0.1,["InnerCircle"],["OuterCircle"],0,0,1,0,0,0.5);