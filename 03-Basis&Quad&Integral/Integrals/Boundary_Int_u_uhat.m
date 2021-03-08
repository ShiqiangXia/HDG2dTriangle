function Buuhat3 = Boundary_Int_u_uhat(k,GQ1DRef_pts,GQ1DRef_wts)
    
    Nuhat = (k+1);
    Nu = (k+2)*(k+1)/2;
    NGQ = length(GQ1DRef_pts);

    Buuhat3 = zeros(Nu,Nuhat,3,numeric_t);
    
    % edge 1:  r=[-1,1], s= -1
    
    % edge 2:  r=[-1,1], s= -r
    
    % edge 3:  r=-1, s= [-1,1]
end