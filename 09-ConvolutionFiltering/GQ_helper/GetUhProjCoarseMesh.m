function uh_coeff = GetUhProjCoarseMesh(k,uh_coarse_GQ_pts,GQ1DRef_pts)
    
    V1D = Vandermonde1D(k,GQ1DRef_pts);
    Lv = (V1D')*V1D;
    Iv = Lv\(V1D');
    IvT = Iv';
    Nk = k+1;
    [~,~,num_ele] = size(uh_coarse_GQ_pts);
    
    uh_coeff = zeros(Nk,Nk,num_ele, numeric_t);
    
    for ii = 1:num_ele
        temp_uh_GQ = uh_coarse_GQ_pts(:,:,ii);
        uh_coeff(:,:,ii) = Iv * temp_uh_GQ * IvT;  
    end
end