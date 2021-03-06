function [M]=ConvolutionFiltering(dom_type,deg_spline,uh,Nx,Ny,mask,Conv_Matrix)
    %
% output M is a 3D array: num_pt x num_pt x num_elemnt
% each row of M(:,:, ele) is the same y different x

% note that this M is the transpose of HDG2D code I wrote before.

%
% uh is the numerical solution (coeffs) on each element
% uh is a matrix  uh_ij = P_i(x)*P_j(y)
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
if  strcmp(dom_type,'Rec')
    
    M = zeros(num_pt,num_pt,Nx*Ny,numeric_t);
    %scale_factor = numeric_t('2')/h;
    
    for m1 = 1:1:Ny % y-direction index
        for  m2 = 1:1:Nx % x-direction index
            if mask(m1,m2) == 0
                % inner domain, we do post-processing
                
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
                M(:,:,(m1-1)*Nx+m2) = temp_sum' ; 
                % each row of M is the same y different x
                % make sure M is consistent with the data structure we use

            end
            
        end
    end
    
elseif strcmp(dom_type,'Lshape')
        
      error('not implemented yet')
      %% L shape

    
    
    
    
end 


end