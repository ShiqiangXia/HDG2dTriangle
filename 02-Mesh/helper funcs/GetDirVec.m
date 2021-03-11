function dir_vec = GetDirVec(Nuhat)
    dir_vec = zeros(Nuhat,1,numeric_t);
    for ii = 0:Nuhat-1
        dir_vec(ii+1) = (-1)^ii;
    end
end