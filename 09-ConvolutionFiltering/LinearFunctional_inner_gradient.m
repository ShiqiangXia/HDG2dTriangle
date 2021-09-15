function [Jh_star, ACh_star,Jh_star_elewise, ACh_star_elewise] = LinearFunctional_inner_gradient...
        ( uH_GQ, grad_uH_x, grad_uH_y, source_f_GQ_pts,...
         vH_GQ, grad_vH_x, grad_vH_y, source_g_GQ_pts,...
         Nx,Ny,mask,GQ1DRef_wts, hx, hy)
     
     % compute two parts
     % part 1: Jh = (uh, g)
     % part 2: ACh = (f,vh)- (grad uh, grad vh)
     
     
     scale_factor = hx * hy * numeric_t('1/4.0');
     
     Jh_star = 0;
     ACh_star = 0;
     
     Jh_star_elewise = zeros(Ny*Nx,1,numeric_t);
     ACh_star_elewise = zeros(Ny*Nx,1,numeric_t);
     
     for ii = 1:Ny
         for jj = 1: Nx
             if mask(ii,jj) == 0
                 temp_ele_id = (ii-1)*Nx+jj;
                 temp_uh_GQ = uH_GQ(:,:,temp_ele_id);
                 temp_g_GQ = source_g_GQ_pts(:,:,temp_ele_id);
                 
                 % (uh, g)
                 func_int = scale_factor*GQ1DRef_wts'*(temp_uh_GQ.*temp_g_GQ)*GQ1DRef_wts;
                 Jh_star_elewise(temp_ele_id,1) = func_int;
                 Jh_star =  Jh_star +func_int;
                     
                 temp_vh_GQ = vH_GQ(:,:,temp_ele_id);
                 temp_f_GQ = source_f_GQ_pts(:,:,temp_ele_id);
                
                 temp_grad_uh_x = grad_uH_x(:,:,temp_ele_id);
                 temp_grad_uh_y = grad_uH_y(:,:,temp_ele_id);
                 
                 temp_grad_vh_x = grad_vH_x(:,:,temp_ele_id);
                 temp_grad_vh_y = grad_vH_y(:,:,temp_ele_id);
                 
                 formula_grad_grad = temp_grad_uh_x.*temp_grad_vh_x ...
                     + temp_grad_uh_y.*temp_grad_vh_y ;
                 
                 % (f,vh) - (grad uh, grad vh)
                 ACh_integral = scale_factor*GQ1DRef_wts'*(temp_vh_GQ.*temp_f_GQ)*GQ1DRef_wts...
                     - scale_factor*GQ1DRef_wts'*(formula_grad_grad)*GQ1DRef_wts;
                 
                 ACh_star_elewise(temp_ele_id,1) = ACh_integral;
                 ACh_star = ACh_star+ACh_integral;
                     
             end
         end
     end
         
end
     
     
     
