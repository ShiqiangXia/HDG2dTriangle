function main(pde_type,varargin)
    
    % example:  
    % main('elliptic','order',1, 'h0',0.2, 'Niter',3, 'refine_flag', 0);
    
    global class_t;
    
    if strcmp(pde_type,'elliptic')
        
        para = SetParameters_Ellipitc(varargin{:});
        
        class_t=para.precision;
        
        EllipticProblemDriver(para);
    elseif strcmp(pde_type,'maxwell')
        
        para = SetParameters_Maxwell(varargin{:});
        
        class_t=para.precision;
        
        MaxwellProblemDriver(para);
        
    else
        error("Wrong PDE type!");
    end
    
end