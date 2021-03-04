function fd = GetDistFunctions(dom_type,x1,y1,x2,y2,bdry_flag)
    
   if strcmp(dom_type,'Rec')
       
       walls = ["bottom","top","left","right"];    
       if isempty(bdry_flag)
           fd = @(p)1;
       else
       
           if ismember( walls(1),bdry_flag)
               f1 = @(p) abs(-p(:,2)+y1);
           else
               f1 = @(p)1;
           end

           if ismember( walls(2),bdry_flag)
               f2 = @(p) abs(p(:,2)-y2);
           else
               f2 = @(p)1;
           end

           if ismember( walls(3),bdry_flag)
               f3 = @(p) abs(x1-p(:,1));
           else
               f3 = @(p)1;
           end

           if ismember( walls(4),bdry_flag)
               f4 = @(p) abs(p(:,1)-x2);
           else
               f4 = @(p)1;
           end

           fd = @(p)min(min(min(f1(p),f2(p)),f3(p)),f4(p));
            
       end
       
   elseif strcmp(dom_type,'L')
       walls = ["bottom","left","right_low","right_high","top_low","top_high"]; 
       
       if isempty(bdry_flag)
           fd = @(p)1;
       else
           if ismember( walls(1),bdry_flag)
               f1 = @(p) abs(p(:,2)-0);
           else
               f1 = @(p)1;
           end

           if ismember( walls(2),bdry_flag)
               f2 = @(p) abs(p(:,1)-0);
           else
               f2 = @(p)1;
           end

           if ismember( walls(3),bdry_flag)
               f3 = @(p) abs(p(:,1)-2);
           else
               f3 = @(p)1;
           end

           if ismember( walls(4),bdry_flag)
               f4 = @(p) abs(p(:,1)-1);
           else
               f4 = @(p)1;
           end
           
           if ismember( walls(5),bdry_flag)
               f5 = @(p) abs(p(:,2)-1);
           else
               f5 = @(p)1;
           end

           if ismember( walls(6),bdry_flag)
               f6 = @(p) abs(p(:,2)-2);
           else
               f6 = @(p)1;
           end
           
           fd = @(p)min ( min(min(min(f1(p),f2(p)),f3(p)),f4(p)),min(f5(p),f6(p)) );
           
       end
       
       
       
           
       
       
       
   end
end