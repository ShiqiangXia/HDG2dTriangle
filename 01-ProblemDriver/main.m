function main(order,h0,Niter)
    
    para = SetParameters(order,h0,Niter);
    
    global class_t;
    class_t=para.precision;
    
    ProblemDriver(para);
end