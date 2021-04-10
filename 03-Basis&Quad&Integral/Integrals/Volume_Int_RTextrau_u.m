function RTuu = Volume_Int_RTextrau_u(k,GQ1DRef_pts,GQ1DRef_wts)
    
    NRTextra = k +1;
    
    Nk_1 = (k+1)*k/2;
    
    RTuu = zeros(NRTextra,Nk_1,2,numeric_t);
    
    NGQ = length(GQ1DRef_pts);
    
    [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts);
    
    V2D = Vandermonde2D(k-1,a_list,b_list);
    
    [V2D_RT1,V2D_RT2] = RTExtraVandermonde2D(k, a_list, b_list);
    
    for ii = 1:NRTextra
        temp_RT1 = V2D_RT1(:,ii); temp_RT1=reshape(temp_RT1,[],NGQ);
        temp_RT2 = V2D_RT2(:,ii); temp_RT2=reshape(temp_RT2,[],NGQ);
        sk = 1;
        for mm= 0:k-1
            for nn = 0:k-1-mm
                temp_u = V2D(:,sk);temp_u=reshape(temp_u,[],NGQ);
                
                RTuu(ii,sk,1) = GQ1DRef_wts'*(temp_RT1.*temp_u.*Jacobian_rs_to_ab)*GQ1DRef_wts;
                RTuu(ii,sk,2) = GQ1DRef_wts'*(temp_RT2.*temp_u.*Jacobian_rs_to_ab)*GQ1DRef_wts;
                
                
                sk = sk+1;
                
            end
        end
    end
    
    
    
end