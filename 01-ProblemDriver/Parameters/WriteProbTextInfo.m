function text_info = WriteProbTextInfo(pb_type,primal_data,adjoint_data)
    pb_type = num2str(pb_type);
    
    if pb_type(2) == 0
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
    elseif pb_type(2) == 1
        text_info = 'eigenvalue problem';
    else
        text_info = ' ';
    end
    
end