
% rectangle mesh
%mymesh1 = Build2DMesh(1,'Rec',0.2,["bottom","top","left"],["right"],0,0,1,1,0);

% L shpaed mesh
mymesh2 = Build2DMesh(0,'L',0.2,["bottom","top_high","right_low","left"],["right_high","top_low"],0);


% Cicle mesh

%mymesh3 = Build2DMesh(0,'Cir',0.05,[],["Circle"],0,0,1);

% Cicle with hole
%mymesh4 = Build2DMesh(0,'CirHole',0.2,["InnerCircle"],["OuterCircle"],0,0,1,0,0,0.5);