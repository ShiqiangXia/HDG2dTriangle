function uh_GQpts = GetUhProjGQpts(uh,k,GQ1DRef_pts)
    % GOAL: Evaluate uh on GQ ref points
    % here uh is defined on square domain with tensor product Qk
    % Output: NGQ x NGQ x num_elements
    % each row: same y different x
    
    NGQ = length(GQ1DRef_pts);
    V1D = Vandermonde1D(k,GQ1DRef_pts);
    [~,~,num_ele] = size(uh);
    
    uh_GQpts = zeros(NGQ, NGQ, num_ele, numeric_t);
    
    
    for ii = 1:num_ele
        uh_GQpts(:,:,ii) = V1D * (uh(:,:,ii)') * (V1D');
        % make sure each row is the same y, different x
    end
    
    
    
end