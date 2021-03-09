function fd = GetDistFunctions(dom_type,bdry_flag,varargin)
    
   if strcmp(dom_type,'Rec')
       
       walls = ["bottom","top","left","right"];
       
       
       if all(ismember(bdry_flag,walls)) ~= 1
           error ("Boundary flags are not correct! Please check the spell.")
           
       end
       
       if isempty(bdry_flag)
           fd = @(p)1;
       else
            x1=varargin{1};
            y1=varargin{2};
            x2=varargin{3};
            y2=varargin{4};

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
       
       if all(ismember(bdry_flag,walls)) ~= 1
           error ("Boundary flags are not correct! Please check the spell.")
           
       end
       if isempty(bdry_flag)
           fd = @(p)1;
       else
           if ismember( walls(1),bdry_flag) % bottom
               f1 = @(p) abs(p(:,2)-0);
           else
               f1 = @(p)1;
           end

           if ismember( walls(2),bdry_flag) % left
               f2 = @(p) abs(p(:,1)-0);
           else
               f2 = @(p)1;
           end

           if ismember( walls(3),bdry_flag) %  right_low
               f3 = @(p) abs(p(:,1)-2);
           else
               f3 = @(p)1;
           end

           if ismember( walls(4),bdry_flag) % right_high
               f4 = @(p) abs(p(:,1)-1) + abs(min(0,p(:,2)-1));
           else
               f4 = @(p)1;
           end
           
           if ismember( walls(5),bdry_flag) % top_low
               f5 = @(p) abs(p(:,2)-1) + abs( min(0, p(:,1)-1) );
           else
               f5 = @(p)1;
           end

           if ismember( walls(6),bdry_flag) % top_high
               f6 = @(p) abs(p(:,2)-2);
           else
               f6 = @(p)1;
           end
           
           fd = @(p)min ( min(min(min(f1(p),f2(p)),f3(p)),f4(p)),min(f5(p),f6(p)) );
           
       end
       
       
   elseif strcmp(dom_type,'Cir')
       
       x1=varargin{1};
       y1=varargin{2};
       r=varargin{3};
       walls =["Circle"];
       
       if isempty(bdry_flag)
           fd = @(p)1;
       else
           if ismember(walls(1),bdry_flag)
               fd = @(p)abs(dcircle(p,x1,y1,r));
           else
               fd = @(p)1;
           end
       end
       
       
   elseif strcmp(dom_type,'CirHole')

        x1=varargin{1};
        y1=varargin{2};
        r1 =varargin{3};

        x2=varargin{4};
        y2=varargin{5};
        r2 =varargin{6};

        walls =["InnerCircle","OuterCircle"];
        
        if isempty(bdry_flag)
           fd = @(p)1;
        else
            if ismember(walls(1),bdry_flag)
               f1 = @(p)abs(dcircle(p,x2,y2,r2));
            else
               f1 = @(p)1;
            end
           
            if ismember(walls(2),bdry_flag)
               f2 = @(p)abs(dcircle(p,x1,y1,r1));
            else
               f2 = @(p)1;
            end
            
            fd = @(p) min(f1(p),f2(p));
            
            
        end
       
       
       
   else
       error("Not implemented yet")
       
       
   end
end