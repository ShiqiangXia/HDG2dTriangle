function [p,e] = Structured2DMesh(dom_type,h0,varargin)
    
    if strcmp(dom_type,'Rec')
        
        x1=varargin{1};
        y1=varargin{2};
        x2=varargin{3};
        y2=varargin{4};
        if nargin-2>=5
            dir = varargin{5};
        else 
            dir = 1; % defaul dirction is (1,1)
        end
        
        lx = (x2-x1);
        ly = (y2-y1);
        
        Nx = round(abs(lx)/h0);
        Ny = round(abs(ly)/h0);
        
        hx = lx/Nx;
        hy = ly/Ny;
        
        p = zeros((Nx+1)*(Ny+1),2);
        e = zeros(2*(Nx*Ny),3);
        
        for jj = 0:Ny
            for ii = 0:Nx
                p(jj*(Nx+1)+ii+1,1) = ii*hx+x1;
                p(jj*(Nx+1)+ii+1,2) = jj*hy+y1;
            end
        end
        
        for jj = 0:Ny-1
            for ii = 1:Nx
                e(jj*Nx+ii,1) = jj*(Nx+1)+ii;
                e(jj*Nx+ii,2) = jj*(Nx+1)+ii+1;
                
                e((Nx*Ny) + jj*Nx+ii,1) = (jj+1)*(Nx+1)+ii;
                e((Nx*Ny) + jj*Nx+ii,2) = (jj+1)*(Nx+1)+ii+1;
                
                if dir == 1 % dirction is (1,1)
                    e(jj*Nx+ii,3) =(jj+1)*(Nx+1)+ii+1;
                    
                    e((Nx*Ny) + jj*Nx+ii,3) = (jj)*(Nx+1)+ii;
                else %direction is (-1,1)
                    e(jj*Nx+ii,3) =(jj+1)*(Nx+1)+ii;
                    
                    e((Nx*Ny) + jj*Nx+ii,3) = (jj)*(Nx+1)+ii+1;
                end
                
            end
        end
        
    elseif strcmp(dom_type,'L')
        
        % structured mesh on [0,2]^2\[1,1]x[2,2]
        
        N = round(1.0/h0);
        
        h = 1.0/N;
        
        p = zeros((2*N+1)*(N+1)+N*(N+1),2);
        e = zeros( 2*(2*N*N+N*N),2);
        
        if nargin-2>=1
            dir = varargin{1};
        else 
            dir = 1; % defaul dirction is (1,1)
        end
        
        % add vertices
        for jj = 0:N
            for ii = 0:2*N
                p(jj*(2*N+1)+ii+1,1) = ii*h;
                p(jj*(2*N+1)+ii+1,2) = jj*h;
            end
        end
        
        for jj = 0:N-1
            for ii = 0:N
                p((2*N+1)*(N+1) + jj*(N+1)+ii+1,1) = ii*h;
                p((2*N+1)*(N+1) + jj*(N+1)+ii+1,2) = (N+jj+1)*h;
            end
        end
        
        
        % add elements
        
        for jj = 0:N-1
            for ii = 1:2*N
                e(jj*2*N+ii,1) = jj*(2*N+1)+ii; 
                e(jj*2*N+ii,2) = jj*(2*N+1)+ii+1; 
                
                e((3*N*N)+jj*2*N+ii,1) = (jj+1)*(2*N+1)+ii; 
                e((3*N*N)+jj*2*N+ii,2) = (jj+1)*(2*N+1)+ii+1;
                
                if dir == 1
                    e(jj*2*N+ii,3) = (jj+1)*(2*N+1)+ii+1;
                    
                    e((3*N*N)+jj*2*N+ii,3) = (jj)*(2*N+1)+ii;
                else
                    e(jj*2*N+ii,3) = (jj+1)*(2*N+1)+ii;
                    
                    e((3*N*N)+jj*2*N+ii,3) = (jj)*(2*N+1)+ii+1;
                end
            end
        end
        
        for jj = 0:N-1
            for ii = 1:N
                if jj == 0
                    
                    e(2*N*N+jj*N+ii,1) = (2*N+1)*N+ii;
                    e(2*N*N+jj*N+ii,2) = (2*N+1)*N+ii+1;
                    
                    e((3*N*N) + 2*N*N+jj*N+ii,1) = (2*N+1)*(N+1)+ii;
                    e((3*N*N) + 2*N*N+jj*N+ii,2) = (2*N+1)*(N+1)+ii+1;
                    
                    if dir == 1
                        e(2*N*N+jj*N+ii,3) = (2*N+1)*(N+1)+ii+1;
                        e((3*N*N) + 2*N*N+jj*N+ii,3) = (2*N+1)*N+ii;
                    else
                        e(2*N*N+jj*N+ii,3) = (2*N+1)*(N+1)+ii;
                        e((3*N*N) + 2*N*N+jj*N+ii,3) = (2*N+1)*N+ii+1;
                    end
                    
                else
                    e(2*N*N+jj*N+ii,1) = (2*N+1)*(N+1) + (jj-1)*(N+1)+ii;
                    e(2*N*N+jj*N+ii,2) = (2*N+1)*(N+1) + (jj-1)*(N+1)+ii+1;
                    
                    e((3*N*N)+ 2*N*N+jj*N+ii,1) = (2*N+1)*(N+1) + (jj)*(N+1)+ii;
                    e((3*N*N)+ 2*N*N+jj*N+ii,2) = (2*N+1)*(N+1) + (jj)*(N+1)+ii+1;
                    
                    if dir == 1 
                        e(2*N*N+jj*N+ii,3) = (2*N+1)*(N+1) + (jj)*(N+1)+ii+1;
                        
                        e((3*N*N)+ 2*N*N+jj*N+ii,3) = (2*N+1)*(N+1) + (jj-1)*(N+1)+ii;
                    else
                        e(2*N*N+jj*N+ii,3) = (2*N+1)*(N+1) + (jj)*(N+1)+ii;
                        
                        e((3*N*N)+ 2*N*N+jj*N+ii,3) = (2*N+1)*(N+1) + (jj-1)*(N+1)+ii+1;
                    end
                    
                end
            end
        end
        
        
        
    else
        error("This type of structured mesh is not implemented yet")
    end
end