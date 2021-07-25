function uh_coeff = GetUhProjCoarseMesh(k,uh_coarse_GQ_pts,GQ1DRef_pts)
    
    % here uh_coeff is a matrix
    % u_{ij} = phi_i(x) * phi_j(y)
    V1D = Vandermonde1D(k,GQ1DRef_pts);
    Lv = (V1D')*V1D;
    Iv = Lv\(V1D');
    IvT = Iv';
    Nk = k+1;
    [~,~,num_ele] = size(uh_coarse_GQ_pts);
    
    uh_coeff = zeros(Nk,Nk,num_ele, numeric_t);
    
    for ii = 1:num_ele
        temp_uh_GQ = uh_coarse_GQ_pts(:,:,ii);
        % make it each row same x diff y, match with uh_coeff data structure
        temp_uh_GQ = temp_uh_GQ'; 
        
        
        uh_coeff(:,:,ii) = Iv * temp_uh_GQ * IvT;  
    end
end