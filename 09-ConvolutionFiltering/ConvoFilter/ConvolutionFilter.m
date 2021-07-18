function ConvolutionFilter(dom_type,k,uh_GQpts,Conv_Matrix, Nx,Ny,N_bd,GQ1DRef_wts)
    % return M = NGQ x NGQ x num_element
    % k is the degree of uh
    % uh_GQpts are uh at GQ pts for each element
    % Nx, Ny are how many elements in each row and col
    [NGQ,~,~] = size(uh_GQpts);
    
    %% so far we assume perioidc bdry so don't consider bound
    if(k==0)
        support =1;  % when use k=0, we need support =1
    else
        support = 2*k; % k>=1 support is in 2*k element
    end
    
%     if N_bd < support
%         error("Not enough neighbors to post-process!\n N_bd should greater than or equal to 2k")
%     end
%     if (N - 2*N_bd) <= 0 
%         error("Not enough inner elements!\n 2*N_bd >= N")
%     end
    
    %% 
    
    
    if strcmp(dom_type,'Rec')
        M = zeros(NGQ,NGQ,Nx*Ny, numeric_t); % Nx*Ny = num_ele
        % go through each element and do convolution
        for m1 = N_bd+1:1:Ny-N_bd % y-direction index
            for  m2 = N_bd+1:1:Nx-N_bd % x-direction index
                sum = zeros(NGQ,NGQ,numeric_t);
                % scan neighbor element with support 2k
                for l1 = m1 - support:1:m1+support % y-direction index
                    for l2 = m2 - support:1:m2+support % x-direction index

                        py = l1-m1;
                        px = l2-m2;
                        
                        idy = mod(l1-1, Ny)+1; % assume periodic
                        idx = mod(l2-1, Nx)+1; % assume periodic

                        scan_ele_id = (idy-1)*Nx+idx ;
                        temp_uh_GQpts = uh_GQpts(:,:,scan_ele_id);
                        temp_uh_GQpts = (temp_uh_GQpts');

                        sum = sum + Conv_Matrix(:,:,px+support+1)...
                            * (temp_uh_GQpts.*GQ1DRef_wts) ...
                            * (Conv_Matrix(:,:,py+support+1)'.*GQ1DRef_wts);

                    end
                end
                
                M(:,:,(m1-1)*Nx+m2) = sum;  
                
            end
        end
    end
end