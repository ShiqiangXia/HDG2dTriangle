function [M]=ConvolutionFiltering(mesh_type,deg_spline,Nk,uh,hx,Nx,hy,Ny,N_bd,Conv_Matrix)
% uh is the numerical solution (coeffs) on each element
% deg_spline : kernel order
% Nk: dimesion of HDG method (polynomial degree Nk-1 for HDG)
% Amtrix, cr: convolution kernel matrix at pts.
% N0: how many elements in a row
% num_bd_elemts: how many elements away from bdary.
%%% pts: quadrature points

[num_pt,~,~] = size(Conv_Matrix);


%% convolution
if(deg_spline==0)
    support =1;  % when use k=0, we need support =1
else
    support = 2*deg_spline; % k>=1 support is in 2*k element
end

% if N_bd < support
%     error("Not enough neighbors to post-process!\n N_bd should greater than or equal to 2k")
% end
% if (N - 2*N_bd) <= 0 
%     error("Not enough inner elements!\n 2*N_bd >= N")
% end

%h = numeric_t('1.0')/N;
if  strcmp(mesh_type,'Rec')
    
    M = zeros(num_pt,num_pt,Nx*Ny,numeric_t);
    %scale_factor = numeric_t('2')/h;
    
    for m1 = N_bd+1:1:Ny-N_bd % y-direction index
        for  m2 = N_bd+1:1:Nx-N_bd % x-direction index
            
            temp_sum = zeros(num_pt,num_pt,numeric_t);
            
            % scan neighbor element with support 2k
            for l1 = m1 - support:1:m1+support % y-direction index
                for l2 = m2 - support:1:m2+support % x-direction index
                    
                    py = l1-m1;
                    px = l2-m2;
                    
                    idy = mod(l1-1, Ny)+1; % assume periodic
                    idx = mod(l2-1, Nx)+1; % assume periodic

                    scan_ele_id = (idy-1)*Nx+idx ;
                    temp_uh = uh(:,:,scan_ele_id);
                    
%                     scan_ele_id = (l1-1)*N+l2 ;
%                     
%                     temp_uh = reshape(uh(:,scan_ele_id),Nk,Nk);
%                     temp_uh = temp_uh';
                
                    temp_sum = temp_sum + Conv_Matrix(:,:,px+support+1) * temp_uh * (Conv_Matrix(:,:,py+support+1)');
                
                     
                end
            end
            M(:,:,(m1-1)*Nx+m2) = temp_sum ; %scale_factor*temp_sum;  
        end
    end
    
elseif strcmp(mesh_type,'Lshape')
        
      error('not implemented yet')
      %% L shape
%     Nin_1 = 2*N - 2*N_bd;
%     Nin_2 = N - 2*N_bd;
      N_elements = 3*N*N;
%     N_inner_elements = Nin_1*Nin_2 + Nin_2*N;
    
    M = zeros(num_pt,num_pt,N_elements,numeric_t);
    scale_factor = numeric_t('2')/h;
    
    for m1 = N_bd+1:1:2*N-N_bd % y-direction index
        for  m2 = N_bd+1:1:2*N-N_bd % x-direction index
            
            if m1>N-N_bd && m2 > N-N_bd
                continue   
            end
            
            temp_sum = zeros(num_pt,num_pt,numeric_t);
            
            % scan neighbor element with support 2k
            for l1 = m1 - support:1:m1+support % y-direction index
                for l2 = m2 - support:1:m2+support % x-direction index
                    
                    py = l1-m1;
                    px = l2-m2;
                    
                    [~,scan_ele_id] = CheckElementIndex(mesh_type,l1,l2,N);
                    
                    temp_uh = reshape(uh(:,scan_ele_id),Nk,Nk);
                    temp_uh = temp_uh';
                
                    temp_sum = temp_sum + Conv_Matrix(:,:,px+support+1) * temp_uh * (Conv_Matrix(:,:,py+support+1)');
                 
                end
            end
            [~,ele_id] = CheckElementIndex(mesh_type,m1,m2,N);
            M(:,:,ele_id) = scale_factor*temp_sum;  
        end
    end
    
    
    
    
end 


end