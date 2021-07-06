function [RT1_ur, RT1_us, RT2_ur, RT2_us] = Volume_Int_RT_du(k_RT, k_u, GQ1DRef_pts,GQ1DRef_wts)
    % Integral of RT basis and gradient of polynomial basis
    % compute (RT_k, du) on ref element
    % k_RT is the degree of RT
    % k_u is the degree of u 
    % There are 4 components we need to compute
    % (RT_1, ur); ?RT_1, us?
    % (RT_2, ur) ; (RT_2, us)
    
    NRT = (k_RT + 3) * (k_RT + 1);
    Nu = (k_u + 1) * (k_u + 2) / 2;
    NGQ = length(GQ1DRef_pts);
    
    RT1_ur = zeros(NRT, Nu, numeric_t);
    RT1_us = zeros(NRT, Nu, numeric_t);
    RT2_ur = zeros(NRT, Nu, numeric_t);
    RT2_us = zeros(NRT, Nu, numeric_t);
    
    [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts);
    
    [V2D_RT1,V2D_RT2] = RTVandermonde2D(k_RT,a_list,b_list); % 1st comp. of RT_k
    
    [V2D_ur,V2D_us] = GradVandermonde2D(k_u,a_list,b_list);
    
    
    for ii = 1:NRT
        temp_RT1 = V2D_RT1(:,ii);
        temp_RT1 = reshape(temp_RT1,[], NGQ);
        temp_RT2 = V2D_RT2(:,ii);
        temp_RT2 = reshape(temp_RT2,[], NGQ);
        ct = 1;
        for mm = 0:k_u
            for nn = 0:k_u - mm
                temp_ur = V2D_ur(:, ct);
                temp_ur = reshape(temp_ur, [], NGQ);
                temp_us = V2D_us(:, ct);
                temp_us = reshape(temp_us, [], NGQ);
                
                RT1_ur(ii, ct) = GQ1DRef_wts'*(temp_RT1.*temp_ur.*Jacobian_rs_to_ab)*GQ1DRef_wts;
                RT1_us(ii, ct) = GQ1DRef_wts'*(temp_RT1.*temp_us.*Jacobian_rs_to_ab)*GQ1DRef_wts;
                RT2_ur(ii, ct) = GQ1DRef_wts'*(temp_RT2.*temp_ur.*Jacobian_rs_to_ab)*GQ1DRef_wts;
                RT2_us(ii, ct) = GQ1DRef_wts'*(temp_RT2.*temp_us.*Jacobian_rs_to_ab)*GQ1DRef_wts;
                
                ct = ct+1;
            end
        end
        
    end
end