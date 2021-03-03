function e = Countclockwise(p,e)
    
    vx = p(:,1);
    vy = p(:,2);
    
    ax = vx(e(:,1));
    ay = vy(e(:,1));
    
    bx = vx(e(:,2));
    by = vy(e(:,2));
    
    cx = vx(e(:,3));
    cy = vy(e(:,3));
    
    % compute the determint of the matrix
    % (ax,ay,1; bx,by,1,cx,cy,1)
    
    D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
    
    i = find(D<0); % <0 clockwise
    
    e(i,:) = e(i,[1 3 2]); % change the ordering 
    
end