function [Auur,Auus] = Volume_Int_u_du(k,GQ1DRef_pts,GQ1DRef_wts)
    
    Nu = (k+2)*(k+1)/2;
    Auur = zeros(Nu,Nu,numeric_t);
    Auus = zeros(Nu,Nu,numeric_t);
    
    % step 1 get all the quad points in 2D
    % (pi,pj)
    NGQ = length(GQ1DRef_pts);
    % Get Gauss Quadpoints on the square
    [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts);
    
    % step 2 get Vandermonde and GradVandermonde matrices.
    
    V2D = Vandermonde2D(k,a_list,b_list);
    
    [V2Dr,V2Ds] = GradVandermonde2D(k,a_list,b_list);
    
    % step 3 do the quadrature.
    
    ct = 1;
    for ii = 0:k
        for jj = 0:k-ii
            temp1 = V2D(:,ct);
            temp1 = reshape(temp1,[],NGQ);
            sk = 1;
            for mm = 0:k
                for nn = 0:k-mm
                    
                    temp2r = V2Dr(:,sk);
                    temp2r = reshape(temp2r,[],NGQ);
                    temp2s = V2Ds(:,sk);
                    temp2s = reshape(temp2s,[],NGQ);
                    
                    Auur(ct,sk) = GQ1DRef_wts'*(temp1.*temp2r.*Jacobian_rs_to_ab)*GQ1DRef_wts;
                    Auus(ct,sk) = GQ1DRef_wts'*(temp1.*temp2s.*Jacobian_rs_to_ab)*GQ1DRef_wts;
                    sk = sk+1;
                end
            end
            
            
            ct = ct+1;
        end
    end
    
    
    
end
