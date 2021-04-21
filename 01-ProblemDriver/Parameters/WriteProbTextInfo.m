function text_info = WriteProbTextInfo(primal_data,adjoint_data)
    
    if primal_data == 0
        text_p = 'u = sin';
    elseif primal_data == 1
        text_p = 'u = corner singular';
    else
        text_p = 'u unknown';
    end
    
    if adjoint_data == 0
        text_a = 'g = sin';
    elseif adjoint_data == 1
        text_a = 'g = 1';
    elseif adjoint_data == 2
        text_a = 'psi discontinuous on boundary';
    else
        text_a ='adjoint data unknown';
    end
    
    text_info = append(text_p,'  ',text_a);
    
end