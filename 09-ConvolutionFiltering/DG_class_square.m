classdef DG_class_square
    
    % a DG solution class for square elements on unit square domain
    % this is used to store uh* obtained from postprocessing
    
    properties
        k % polynomial degree
        %mesh % the mesh where DG solution is defined
        hx 
        hy
        Nx
        Ny
        GQ1DRef_pts
        mask % this is used to restrict DG solution on some elements
        dg_gq_pts % DG solution at GQ points for each element
        dg_coeff % DG soluttion coefficient
        
    end
    methods
        function obj = DG_class_square(k,hx,hy,Nx,Ny,GQ1DRef_pts,mask)
            obj.k = k;
            %obj.mesh = mesh;
            obj.hx = hx;
            obj.hy = hy;
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.mask = mask;
            obj.GQ1DRef_pts = GQ1DRef_pts;
        end
        
        function obj = set_dg_coeff(dg_coeff)
            obj.dg_coeff = dg_coeff;
        end
        
        function obj = set_dg_gq_pts(obj,dg_gq_pts)
            obj.dg_gq_pts = dg_gq_pts;
            obj = obj.get_coeff_from_gq();
        end
        
        function obj=get_coeff_from_gq(obj)
            %  V1D * coeff' * V1D' = GQ_mat
            V1D = Vandermonde1D(obj.k,obj.GQ1DRef_pts);
            Lv = (V1D') * V1D;
            Iv = Lv\(V1D');
            IvT = Iv';
            Nk = obj.k +1;
            Nele = obj.Nx * obj.Ny;
            coeff = zeros(Nk,Nk,Nele,numeric_t);
            for i = 1:obj.Ny
                for j = 1:obj.Nx
                    if obj.mask(i,j) == 0
                        temp_mat = obj.dg_gq_pts(:,:,(i-1)*obj.Nx+j);
                        
                        coeff(:,:,(i-1)*obj.Nx+j) = (Iv * temp_mat * IvT)';
                        
                    end
                end
            end
            
            obj.dg_coeff = coeff;
        end
        
        function res = eval(obj,xpts, ypts)
            % evaluate DG solution at points (x,y) in (xpts, ypts)
            n = length(xpts);
            if n ~= length(ypts)
                error('Input xpts and ypts must have the same length!')
            end
            
            res = zeros(n,1,numeric_t);
            tol = 1e-8;
           
            for i = 1:n
                xpt = xpts(i);
                ypt = ypts(i);
                if xpt>1+tol || xpt<0-tol || ypt>1+tol || ypt<0-tol
                    print('Invalid point! Points should be in the unit square')
                end
                % step 1: determine which element xpt and ypt belongs to
                if mod(xpt, obj.hx)<tol
                    % point is on the board face
                    temp = round(xpt/obj.hx);
                    x_id = [temp, temp+1];
                else
                    x_id = floor(xpt/obj.hx)+1;
                end
                
                if mod(ypt, obj.hy)<tol
                    % point is on the board face
                    temp = round(ypt/obj.hy);
                    y_id = [temp, temp+1];
                else
                    y_id = floor(ypt/obj.hy)+1;
                end
                
                
                % step 2: transfer to ref point
                 
                % step 3: evaluate
                ct = 0;
                temp_val = 0;
                for s = 1:length(y_id)
                    for t = 1:length(x_id)
                        xind = x_id(t);
                        yind = y_id(s);
                        if xind>0 && xind <= obj.Nx && yind >0 && yind <= obj.Ny
                            if obj.mask(yind,xind) == 0
                                ct = ct + 1;
                                xref = (xpt - (xind-0.5)*obj.hx)/(0.5*obj.hx);
                                yref = (ypt - (yind-0.5)*obj.hy)/(0.5*obj.hy);
                                V1D_x = Vandermonde1D(obj.k,xref);
                                V1D_y = Vandermonde1D(obj.k,yref);                                
                                temp_val = temp_val + V1D_x * obj.dg_coeff(:,:,(yind-1)*obj.Nx+xind) * V1D_y';
                            end
                        end
                    end
                end
                if ct ~= 0
                    temp_val = temp_val / ct; % average for points on the faces
                end
                
                res(i,1) = temp_val;
                
                
                
                
               
            end
        end
        
        function [grad_uh_x, grad_uh_y] = get_DG_grad_pts(obj, x_pts, y_pts)
            % Get the grad uh at points (x_pts, y_pts)
            if length(x_pts) ~= length(y_pts)
                error('x_pts and y_pth lengths dont match!')
            end
            
            poly_k = obj.k;
            Lv_x = Vandermonde1D(poly_k,x_pts);
            Lv_y = Vandermonde1D(poly_k,y_pts);
            Dv_x = GradVandermonde1D(poly_k,x_pts);
            Dv_y = GradVandermonde1D(poly_k,y_pts);
            Npt = length(x_pts);
            [~,~,Nele] = size(obj.dg_coeff);
            grad_uh_x = zeros(Npt,Npt,Nele);
            grad_uh_y = zeros(Npt,Npt,Nele);
            for ii = 1:obj.Ny
                for jj = 1: obj.Nx
                    if obj.mask(ii,jj) == 0
                        temp_ele_id = (ii-1)*obj.Nx+jj;
                        temp_uh_coeff = obj.dg_coeff(:,:,temp_ele_id);
                        temp_grad_uh_x = Dv_x * temp_uh_coeff * Lv_y' * 2.0/obj.hx;
                        temp_grad_uh_y = Lv_x * temp_uh_coeff * Dv_y'* 2.0/obj.hy;
                        grad_uh_x(:,:,temp_ele_id) = temp_grad_uh_x';
                        grad_uh_y(:,:,temp_ele_id) = temp_grad_uh_y';

                    end
                end
            end
            
        end
    end
end