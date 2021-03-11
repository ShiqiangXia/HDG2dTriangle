function ord = GetOrder(dof,err)
    
    nn = size(dof,1);
    ord = zeros(nn,1);
    
    top = log(err(1:end-1,1)./err(2:end,1));
    bot = 0.5*log(dof(2:end,1)./dof(1:end-1,1));
    
    
    ord(2:end,1) = top./bot;
    
end