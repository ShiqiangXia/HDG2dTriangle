function Buuhat3 = Boundary_Int_u_uhat(k,GQ1DRef_pts,GQ1DRef_wts)
    
    % Notice that here uhat has the same orientation as the triangle(counterclock wise)
    Nuhat = (k+1);
    Nu = (k+2)*(k+1)/2;
    NGQ = length(GQ1DRef_pts);

    Buuhat3 = zeros(Nu,Nuhat,3,numeric_t);
    
    V1D = Vandermonde1D(k,GQ1DRef_pts);% legendre polynomail at GQ points
    
    % edge 1:  r=[-1,1], s= -1 , V1---> V2
    
    r1 = GQ1DRef_pts;
    s1 = -1*ones(NGQ,1,numeric_t);
    
    [a1,b1] = RStoAB(r1,s1);
    
    V2D1 = Vandermonde2D(k,a1,b1);
    
    ct = 1;
    for ii = 0:k
        for jj = 0:k-ii
            
            temp1 = V2D1(:,ct);
            
            for nn = 1:k+1
                temp2 = V1D(:,nn);
                Buuhat3(ct,nn,1) = GQ1DRef_wts'*(temp1.*temp2);
            end
            
            ct = ct+1;
        end
    end
    
    
    
    % edge 2:  r=-s, s=[-1,1]  V2--->V3
    
    s2 = GQ1DRef_pts;
    r2 = -s2;
    
    [a2,b2] = RStoAB(r2,s2);
    
    V2D2 = Vandermonde2D(k,a2,b2);
    %scale_factor = numeric_t('sqrt(2)');
    ct = 1;
    for ii = 0:k
        for jj = 0:k-ii
            temp1 = V2D2(:,ct);
            for nn = 1:k+1
                temp2 = V1D(:,nn);
                Buuhat3(ct,nn,2) = GQ1DRef_wts'*(temp1.*temp2);
                
            end
            ct = ct+1;
            
        end
    end
    
    % edge 3:  r=-1, s= [1,-1]  V3--->V1
    
    r3 = -1*ones(NGQ,1,numeric_t);
    s3 = - GQ1DRef_pts;
    
    
    [a3,b3] = RStoAB(r3,s3);
    
    V2D3 = Vandermonde2D(k,a3,b3);
    
    ct = 1;
    for ii = 0:k
        for jj = 0:k-ii

            temp1 = V2D3(:,ct);
            for nn = 1:k+1
                temp2 = V1D(:,nn);
                Buuhat3(ct,nn,3) = GQ1DRef_wts'*(temp1.*temp2);
                
            end
            ct = ct+1;
        end
    end
    
    
end