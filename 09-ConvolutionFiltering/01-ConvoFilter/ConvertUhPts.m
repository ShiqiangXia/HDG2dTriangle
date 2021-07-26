function  M1 = ConvertUhPts(M,k,source_pts,target_pts)
    
    V_source = Vandermonde1D(k,source_pts);
    V_target = Vandermonde1D(k,target_pts);
    Npt = length(target_pts);
    Lv = V_target/V_source;
    
    [~,~,num_ele] = size(M);
    
    M1 = zeros(Npt,Npt,num_ele);
    
    for ii = 1:num_ele
        M1(:,:,ii) = Lv * M(:,:,ii) * Lv';
    end
    
end