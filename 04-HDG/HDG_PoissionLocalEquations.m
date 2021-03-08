function [N_1,N_2,N_3,M_Loc] = HDG_PoissionLocalEquations(...
        Jk,vertice_list,tau,Aqq,Auur,Auus,Auu3,Buuhat3)
    
    % Jk is the jacobian of the element
    % vertice_list [x1,y1;x2,y2;x3,y3]--?counterclockwise as the element vertice index
    % HDG tau
    
    %---------------------------------------------------------------------
    % M_loc
    % (q,r)
    Mqq = Jk*Aqq;
    
    %(u,div.r)
    % get the affine map and it's inverse
    
    [~,Inv_AffineMap] = Ref_Tri_Map(Jk,vertice_list);
    
                        
    temp1 = Jk*(Auur*Inv_AffineMap(1,1) + Auus*Inv_AffineMap(1,2));
    temp2 = Jk*(Auur*Inv_AffineMap(2,1) + Auus*Inv_AffineMap(2,2));
    
    Muq = [temp1,temp2];
    
    % <tau u, w>
    [e1,e2,e3,n1,n2,n3] = GetTriFaceInfo(V1,V2,V3);
    
    Muu = 0.5*tau*(e1*Auu3(:,:,1)+e2*Auu3(:,:,2)+e3*Auu3(:,:,3));
    
    M_Loc = [Mqq, -Muq'; Muq,Muu ];
    
    %---------------------------------------------------------------------
    % N_1,N_2,N_3
    
    N_1=[-n1(1)*Buuhat3(:,:,1);...
        -n1(2)*Buuhat3(:,:,1);...
        0.5*tau*e1*Buuhat3(:,:,1)];
    
    N_2=[-n2(1)*Buuhat3(:,:,2);...
        -n2(2)*Buuhat3(:,:,2);...
        0.5*tau*e2*Buuhat3(:,:,2)];
    
    N_3=[-n3(1)*Buuhat3(:,:,3);...
        -n3(2)*Buuhat3(:,:,3);...
        0.5*tau*e3*Buuhat3(:,:,3)];
    
    
end