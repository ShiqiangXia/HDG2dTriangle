function [e1,e2,e3,n1,n2,n3] = GetTriFaceInfo(V1,V2,V3)
    % compute the length of the edeges and the normal vector
    % ASSUME: V1,V2,V3 are labeled counterclockwise
    % Vi = [xi;yi] colomum vector
    % 
    % if not counterclockwise, we just need to check the determint 
    %      | x1,y1,1 |
    % det =| x2,y2,1 |
    %      | x3,y3,1 |
    % if det>0 counterclockwise
    % if det<0 clockwise
    % then n1 = - sign(det)*R*(V2-V1), ect.
    
    e1 = VectorNorm(V2-V1);
    e2 = VectorNorm(V3-V2);
    e3 = VectorNorm(V1-V3);
    
    Rot_mat = [0,-1;1,0]; % counterclock wise rotate 90 degree
    
    n1 = - Rot_mat*(V2-V1);
    n1 = n1/VectorNorm(n1);
    
    n2 = - Rot_mat*(V3-V2);
    n2 = n2/VectorNorm(n2);
    
    n3 = - Rot_mat*(V1-V3);
    n3 = n3/VectorNorm(n3);
    
end