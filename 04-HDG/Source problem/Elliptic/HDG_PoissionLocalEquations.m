function [Ns,M_Loc] = HDG_PoissionLocalEquations(...
        Jk,vertice_list,tau,Auu,Auur,Auus,Auu3,Buuhat3,uhat_dir_list)
    
    % Jk is the jacobian of the element
    % vertice_list [x1,y1;x2,y2;x3,y3]--?counterclockwise as the element vertice index
    % HDG tau
    
    [Nu,Nuhat,~] = size(Buuhat3);
    Nq = 2*Nu;
    Aqq = blkdiag(Auu,Auu);
    
    dir_vec = GetDirVec(Nuhat); % correct the uhat oritation 
    dir_vec = dir_vec';
    % uh_hat has to have a unqiue oritation, we use the oritentiaon:
    % small_verice_idx ---> large_vertice_idx
    % the local matrix Buuhat3 = <uhat,u> is using the counterclock wise
    % orientation. it's not unique in two triangles, so we need to fix it. 
    % the difference is legendre_i(-x) = (-1)^i*legendre_i(x)
    
    %---------------------------------------------------------------------
    % M_loc
    % (q,r)
    Mqq = Jk*Aqq;
    
    %(u,div.r)
    % get the affine map and it's inverse
    
    [~,Inv_AffineMap] = Ref_Tri_Map(Jk,vertice_list);
    
                        
    temp1 = Jk*(Auur*Inv_AffineMap(1,1) + Auus*Inv_AffineMap(2,1));
    temp2 = Jk*(Auur*Inv_AffineMap(1,2) + Auus*Inv_AffineMap(2,2));
    
    Muq = [temp1,temp2];
    
    % <tau u, w>
    [e_list,n1,n2,n3] = GetTriFaceInfo(vertice_list);
    
    Muu = 0.5*tau*(e_list(1)*Auu3(:,:,1)+e_list(2)*Auu3(:,:,2)+e_list(3)*Auu3(:,:,3));
    
    M_Loc = [Mqq, -Muq'; Muq,Muu ];
    
    %---------------------------------------------------------------------
    % N_1,N_2,N_3
    %scale_factor = numeric_t('sqrt(2)'); 
    
    
    
    N_1=0.5*e_list(1)* [-n1(1)*Buuhat3(:,:,1);...
        -n1(2)*Buuhat3(:,:,1);...
        tau*Buuhat3(:,:,1)];
    
    if uhat_dir_list(1,1) == 0  % reverse uhat direction
        
        N_1 = N_1.*dir_vec;
    end
    
    N_2=0.5*e_list(2) *[-n2(1)*Buuhat3(:,:,2);...
        -n2(2)*Buuhat3(:,:,2);...
        tau*Buuhat3(:,:,2)];
    
    if uhat_dir_list(1,2) == 0  % reverse uhat direction
        
        N_2 = N_2.*dir_vec;
    end
    
    N_3= 0.5*e_list(3)*[-n3(1)*Buuhat3(:,:,3);...
        -n3(2)*Buuhat3(:,:,3);...
        tau*Buuhat3(:,:,3)];
    
    if uhat_dir_list(1,3) == 0  % reverse uhat direction
        N_3 = N_3.*dir_vec;
    end
    
    
    Ns = zeros(Nq+Nu,Nuhat,3,numeric_t);
    Ns(:,:,1) = N_1;
    Ns(:,:,2) = N_2;
    Ns(:,:,3) = N_3;
    
end