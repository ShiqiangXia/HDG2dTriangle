function main(pde_type,varargin)
    
    % example:  
    
    % main('maxwell','order',1, 'h0',0.2, 'Niter',3,'refine_flag', 0);
    
    % main('elliptic','order',3, 'h0',0.2, 'Niter',3, 'refine_flag', 0,'pb_type',2011,'dom_type','Rec','primal',0,'adjoint',0,'post_process_flag',1,'err_cal_flag',1)
    
    global class_t;
    
    if strcmp(pde_type,'elliptic')
        
        para = SetParameters_Ellipitc(varargin{:});
        
        class_t=para.precision;
        
        EllipticProblemDriver(para);
        
    elseif strcmp(pde_type,'maxwell')
        
        para = SetParameters_Maxwell(varargin{:});
        
        class_t=para.precision;
        
        MaxwellProblemDriver(para);
        
    elseif strcmp(pde_type,'elliptic_trick')
        
        para = SetParameters_Ellipitc(varargin{:});
        
        class_t=para.precision;
        
        EllipitcProblemDriver_new_trick(para);
        
        
    else
        error("Wrong PDE type!");
    end
    
end