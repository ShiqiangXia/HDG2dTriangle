function ReportProblem(para)
    
    %% Get infor
    pb_type = num2str(para.pb_type);
    if strcmp(pb_type(1),'1')
        pb1 = 'PDE';
    elseif strcmp(pb_type(1),'2')
        pb1 = 'Functional';
    else
        
    end
    
    if strcmp(pb_type(2),'0')
        pb2 = 'source problem';
    elseif strcmp(pb_type(2),'1')
        pb2 = 'eigen problem';
    else   
    end
    
    if strcmp(pb_type(3),'1')
        if strcmp(pb_type(2),'0')
            pb_eq = 'Eq: -Laplace u = f\n   u = uD on bounary\n';
        elseif strcmp(pb_type(2),'1')
            pb_eq = 'Eq: -Laplace u = lam * u\n   u = 0 on bounary\n';
        end
    elseif strcmp(pb_type(3),'5')
        if strcmp(pb_type(2),'0')
            pb_eq = 'Eq: mu*w - curl(u)=0\n    curl(w)- eps*omg^2*u + eps*grad(p) = j\n    div(eps*u)=0\n    u x n = g on boundary\n    p = 0 on boundary\n';
        elseif strcmp(pb_type(2),'1')
            pb_eq = 'Eq: mu*w - curl(u)=0\n    curl(w)- eps*omg^2*u + eps*grad(p) = j\n    div(eps*u)=0\n    u x n = g on boundary\n    p = 0 on boundary\n';
        end
    else  
    end
    
    if strcmp(pb_type(4),'0')
        func = '';
    elseif strcmp(pb_type(4),'1')
        func = 'J(u) = (g,u)';
    elseif strcmp(pb_type(4),'2')
        func = 'J(u) = <q*n,phi>';
    elseif strcmp(pb_type(4),'3')
        func = 'J(u) = eigenvalue(u)';
    else
    end
    
    if para.refine_flag == 0
        ref='Uniform';
    elseif para.refine_flag == -1
        ref ='Remesh with half h';  
    elseif para.refine_flag == 1
        ref='Adaptive with RGB';
    elseif para.refine_flag == 2
        ref='Adaptive with RG';
    elseif para.refine_flag == 3
        ref='Adaptive with NVB';
    else
    end
    
    %% Print results
    fprintf('\n')
    fprintf('---------- Report ------------\n')
    
    fprintf("%s %s\n",pb1,pb2 );
    fprintf(pb_eq);
    if ~strcmp(pb_type(4),'0')
        fprintf("Functional: %s\n",func);
    end
    
    fprintf('Dom type: %s\n', para.dom_type);
    fprintf('order k = %d\n', para.order);
    fprintf('Refinement: %s\n',ref);
    
    
        
    
    
    
end