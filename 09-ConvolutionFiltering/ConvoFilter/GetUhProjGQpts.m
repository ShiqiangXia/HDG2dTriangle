function uh_pro_GQpts = GetUhProjGQpts(uh_coeff,k,GQ1DRef_pts)
    % Output: NGQ x NGQ x num_elements
    % each row: same y different x
    % plot
    NGQ = length(GQ1DRef_pts);
    V1D = Vandermonde1D(k,GQ1DRef_pts);
    [~,~,num_ele] = size(uh_coeff);
    
    uh_pro_GQpts = zeros(NGQ, NGQ, num_ele, numeric_t);
    
    
    for ii = 1:num_ele
        uh_pro_GQpts(:,:,ii) = V1D * (uh_coeff(:,:,ii)') * (V1D');
        % make sure each row is the same y, different x
    end
    
    
    
end