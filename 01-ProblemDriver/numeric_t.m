function r = numeric_t(expression)
    global class_t;
    if (nargin > 0)
        if(strcmpi(class_t,'mp')), r = mp(expression); 
        else
            if isnumeric(expression)
                r = expression;
            else
                r = eval(expression);
            end
        end
    else
        r = class_t;
    end
end