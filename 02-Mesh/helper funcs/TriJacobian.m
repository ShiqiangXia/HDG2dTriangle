function J_list = TriJacobian(p,e)
    ne = size(e,1);
    if isempty(numeric_t)
        J_list = zeros(ne,1);
    else
    J_list  = zeros(ne,1,numeric_t);
    end
    
    for ii = 1:ne
        V1 = p(e(ii,1),:)';
        V2 = p(e(ii,2),:)';
        V3 = p(e(ii,3),:)';
        
        m = numeric_t('0.5')*[V2-V1,V3-V1];
        J_list(ii,1) = det(m);
    end
end