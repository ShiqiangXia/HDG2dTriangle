function [GQ_x, GQ_y] = GetPhyGQPts(structured_flag,dom_type,h0,GQ1DRef_wts, varargin)
    
    if structured_flag ~= 1
        error('None structured type is not implemented')
    end
    
    
    NGQ = length(GQ1DRef_wts);
    
    if strcmp(dom_type,'Rec')
        x1=varargin{1};
        y1=varargin{2};
        x2=varargin{3};
        y2=varargin{4};
        
        lx = (x2-x1);
        ly = (y2-y1);
        
        Nx = round(abs(lx)/h0);
        Ny = round(abs(ly)/h0);
        
        hx = lx/Nx;
        hy = ly/Ny;
        
        GQ_x = zeros(NGQ, Nx*Ny, numeric_t);
        GQ_y = zeros(NGQ, Nx*Ny, numeric_t);
        
        for jj = 0:Ny-1
            for ii = 1:Nx
                xmid = (ii-1)*hx + hx*0.5;
                GQ_x(:, jj*Nx + ii) = hx*0.5 * GQ1DRef_wts + xmid;
                
                ymid = jj*hy + hy*0.5;
                GQ_y(:,jj*Nx+ii) = hy*0.5 * GQ1DRef_wts + ymid;
            end
        end
        
    else
        error('Dom type has not implemented yet')
    end
end