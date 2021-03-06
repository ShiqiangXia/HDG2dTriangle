function Auu3 = Boundary_Int_u_u(k,GQ1DRef_pts,GQ1DRef_wts)
    
    Nu = (k+2)*(k+1)/2;
    
    NGQ = length(GQ1DRef_pts);
    
    Auu3 = zeros(Nu,Nu,3,numeric_t);
    
    
    % edge 1:  r=[-1,1], s= -1
    r1 = GQ1DRef_pts;
    s1 = -1*ones(NGQ,1,numeric_t);
    
    [a1,b1] = RStoAB(r1,s1);
    
    V2D1 = Vandermonde2D(k,a1,b1);
    
    ct = 1;
    for ii = 0:k
        for jj = 0:k-ii
            sk = 1;
            temp1 = V2D1(:,ct);
            for mm = 0:k
                for nn = 0:k-mm
                    temp2 = V2D1(:,sk);
                    
                    Auu3(ct,sk,1) = GQ1DRef_wts'*(temp1.*temp2);
                    
                    sk=sk+1;
                end
            end
            
            ct = ct+1;
        end
    end
    
    
    
    % edge 2:  r=[-1,1], s= -r
    r2 = GQ1DRef_pts;
    s2 = -r2;
    [a2,b2] = RStoAB(r2,s2);
    
    V2D2 = Vandermonde2D(k,a2,b2);
    ct = 1;
    for ii = 0:k
        for jj = 0:k-ii
            sk = 1;
            temp1 = V2D2(:,ct);
            for mm = 0:k
                for nn = 0:k-mm
                    temp2 = V2D2(:,sk);
                    
                    Auu3(ct,sk,2) = GQ1DRef_wts'*(temp1.*temp2);
                    
                    sk=sk+1;
                end
            end
            
            ct = ct+1;
        end
    end
    
    
    % edge 3:  r=-1, s= [-1,1]
    
    r3 = -1*ones(NGQ,1,numeric_t);
    s3 = GQ1DRef_pts;
    
    
    [a3,b3] = RStoAB(r3,s3);
    
    V2D3 = Vandermonde2D(k,a3,b3);
    
    ct = 1;
    for ii = 0:k
        for jj = 0:k-ii
            sk = 1;
            temp1 = V2D3(:,ct);
            for mm = 0:k
                for nn = 0:k-mm
                    temp2 = V2D3(:,sk);
                    
                    Auu3(ct,sk,3) = GQ1DRef_wts'*(temp1.*temp2);
                    
                    sk=sk+1;
                end
            end
            
            ct = ct+1;
        end
    end
    
    
    
end