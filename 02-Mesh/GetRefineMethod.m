function rlt = GetRefineMethod(flag)
    if flag == 1
        rlt = 'RGB';
    elseif flag == 2
        rlt = 'RG';
    elseif flag == 3
        rlt = 'NVB';
    else
        error('Refine method %d is not supported.',flag)
    end
end