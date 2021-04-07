function RTuuhat3 = Boundary_Int_RTExtra_u_uhat(k,GQ1DRef_pts,GQ1DRef_wts)
    
    % compute the extra RT basis <P_{i+1,k-i}, mu> and <P_{i,k+1-i}, mu>
    
    % Notice that here uhat has the same orientation as the triangle(counterclock wise)
    Nuhat = (k+1);
    Nu = (k+2)*(k+1)/2;
    NGQ = length(GQ1DRef_pts);

    RTuuhat3 = zeros(Nuhat,Nuhat,3,2,numeric_t);
    
    V1D = Vandermonde1D(k,GQ1DRef_pts);% legendre polynomail at GQ points
    
    % edge 1:  r=[-1,1], s= -1 , V1---> V2
    
    r1 = GQ1DRef_pts;
    s1 = -1*ones(NGQ,1,numeric_t);
    
    [a1,b1] = RStoAB(r1,s1);
    
    [V2D1_1,V2D1_2] = RTExtraVandermonde2D(k,a1,b1);
    
    
    for ct = 1:k+1
        tempRT1 = V2D1_1(:,ct);
        tempRT2 = V2D1_2(:,ct);
        for nn = 1:k+1
            temp2 = V1D(:,nn);
            RTuuhat3(ct,nn,1,1) = GQ1DRef_wts'*(tempRT1.*temp2);
            RTuuhat3(ct,nn,1,2) = GQ1DRef_wts'*(tempRT2.*temp2);

        end
    end
    
    
    
    % edge 2:  r=-s, s=[-1,1]  V2--->V3
    
    s2 = GQ1DRef_pts;
    r2 = -s2;
    
    [a2,b2] = RStoAB(r2,s2);
    
    [V2D2_1,V2D2_2] = RTExtraVandermonde2D(k,a2,b2);
    
    %scale_factor = numeric_t('sqrt(2)');
   
    for ct = 1:k+1
        
        tempRT1 = V2D2_1(:,ct);
        tempRT2 = V2D2_2(:,ct);
        for nn = 1:k+1
            temp2 = V1D(:,nn);
            RTuuhat3(ct,nn,2,1) = GQ1DRef_wts'*(tempRT1.*temp2);
            RTuuhat3(ct,nn,2,2) = GQ1DRef_wts'*(tempRT2.*temp2);

        end
       
    end
    
    % edge 3:  r=-1, s= [1,-1]  V3--->V1
    
    r3 = -1*ones(NGQ,1,numeric_t);
    s3 = - GQ1DRef_pts;
    
    
    [a3,b3] = RStoAB(r3,s3);
    
    [V2D3_1,V2D3_2] = RTExtraVandermonde2D(k,a3,b3);
    
  
    for ct = 1:k+1
        
        tempRT1 = V2D3_1(:,ct);
        tempRT2 = V2D3_2(:,ct);
        for nn = 1:k+1
            temp2 = V1D(:,nn);
            RTuuhat3(ct,nn,3,1) = GQ1DRef_wts'*(tempRT1.*temp2);
            RTuuhat3(ct,nn,3,2) = GQ1DRef_wts'*(tempRT2.*temp2);

        end
      
    end
    
    
end