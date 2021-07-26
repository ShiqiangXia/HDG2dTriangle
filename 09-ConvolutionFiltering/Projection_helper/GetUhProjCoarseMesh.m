function uh_coeff = GetUhProjCoarseMesh(k,uh_coarse_GQ_pts,GQ1DRef_pts)
    % GOAL: Interpolation to Qk space based on uh at GQ points
    % OUTPUT: uh_coeff is a matrix
    % u_{ij} = phi_i(x) * phi_j(y)
    
    % base on uh at GQ poits to do a least square approximation
    % namely find uh in Qk such that uh 
    % minimize sum (uh(pi) -uh_coarse_GQ(pi))^2
    
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