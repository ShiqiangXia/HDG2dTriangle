function PercentPlot(s,j)
    
    total = sum(s);
    nn = size(s,1);
    mylist = zeros(nn,1);
    for ii = 1:nn
        mylist(ii) = sum(s(1:ii))/total;
        %mylist(ii) = s(ii)/total;
    end
    
    figure;
    plot(1:nn,mylist,'x--');
    xline(j,'--r');
    
    %eq = 1/nn;
    %yline(eq,'--r');
    title('percent vs N');
end