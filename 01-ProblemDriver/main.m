function main()
    
    para = SetParameters();
    
    global class_t;
    class_t=para.precision;
    
    ProblemDriver(para);
end