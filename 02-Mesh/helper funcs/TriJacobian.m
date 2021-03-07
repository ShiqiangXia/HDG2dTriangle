function J_list = TriJacobian(p,e)
    ne = length(e);
    J_list  = zeros(ne,1,numeric_t);
    
    for ii = 1:ne
        V1 = p(e(ii,1),:)';
        V2 = p(e(ii,2),:)';
        V3 = p(e(ii,3),:)';
        
        m = numeric_t('0.5')*[V2-V1,V3-V1];
        J_list(ii) = det(m);
    end
end