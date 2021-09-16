function  err_inner_estimate = FuncErrInnerEstimate(...
                    mask,Nx,Ny,hx,hy,...
                    q1_x,q1_y,...
                    q2_x,q2_y,...
                    q3_x,q3_y,...
                    q4_x, q4_y,...
                    GQ1DRef_wts...
                    )
       scale_factor = hx * hy * numeric_t('1/4.0');
       err_inner_estimate = 0.0 ;      
       for ii = 1:Ny
         for jj = 1: Nx
             if mask(ii,jj) == 0
                 ele_id = (ii-1)*Nx+jj;
                 diff1 = (q1_x(:,:,ele_id) - q2_x(:,:,ele_id)).*(q3_x(:,:,ele_id) - q4_x(:,:,ele_id));
                 diff2 = (q1_y(:,:,ele_id) - q2_y(:,:,ele_id)).*(q3_y(:,:,ele_id) - q4_y(:,:,ele_id));
                 formula = diff1 + diff2;
                 
                 err_inner_estimate = err_inner_estimate +...
                     scale_factor*GQ1DRef_wts'*(formula)*GQ1DRef_wts;
                 
             end
         end
       end
                
end