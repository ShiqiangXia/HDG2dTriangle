function ord = GetOrderH(h_list,err)
    
    nn = size(h_list,1);
    ord = zeros(nn,1);
    
    top = log(err(1:end-1,1)./err(2:end,1));
    bot = -log(h_list(2:end,1)./h_list(1:end-1,1));
    
    
    ord(2:end,1) = top./bot;
    
end