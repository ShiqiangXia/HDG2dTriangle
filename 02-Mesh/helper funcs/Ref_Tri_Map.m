function [AffineMap,Inv_AffineMap]=Ref_Tri_Map(Jk,vertice_list)
    % compute the affine map from the ref element to any
    % triangle difined by the vertice_list
    %[x,y]' = a_mat *[r+1,s+1]'+[x1,y1]'
    % inv_a_mat = Inverse(a_mat)
    
    V1 = vertice_list(1,:)';
    V2 = vertice_list(2,:)';
    V3 = vertice_list(3,:)';
    AffineMap = numeric_t('0.5')*[V2-V1,V3-V1];
    Inv_AffineMap = (1.0/Jk) * [AffineMap(2,2),-AffineMap(1,2);...
                            -AffineMap(2,1),AffineMap(1,1)];
end
