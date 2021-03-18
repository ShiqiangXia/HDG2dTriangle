function main(varargin)
    
    % example:  main('order',1, 'h0',0.2, 'Niter',3, 'refine_flag', 0);
    
    para = SetParameters(varargin{:});
    
    global class_t;
    class_t=para.precision;
    
    ProblemDriver(para);
end