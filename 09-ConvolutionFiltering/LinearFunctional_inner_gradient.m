function [Jh_star, ACh_star,Jh_star_elewise, ACh_star_elewise] = LinearFunctional_inner_gradient...
        ( func_type,pde_ype,...
        uH_GQ, grad_uH_x, grad_uH_y, source_f_GQ_pts,...
         vH_GQ, grad_vH_x, grad_vH_y, source_g_GQ_pts,...
         Nx,Ny,mask,GQ1DRef_wts, hx, hy)
     
     % compute two parts
     % part 1: term_A = (uh, g)
     % part 2: term_B = (f,vh)- (grad uh, grad vh)
     
     
     scale_factor = hx * hy * numeric_t('1/4.0');
     
     term_A = 0;
     term_B = 0;
     
     term_A_elewise = zeros(Ny*Nx,1,numeric_t);
     term_B_elewise = zeros(Ny*Nx,1,numeric_t);
     
     for ii = 1:Ny
         for jj = 1: Nx
             if mask(ii,jj) == 0
                 temp_ele_id = (ii-1)*Nx+jj;
                 temp_uh_GQ = uH_GQ(:,:,temp_ele_id);
                 temp_g_GQ = source_g_GQ_pts(:,:,temp_ele_id);
                 
                 % (uh, g)
                 func_int = scale_factor*GQ1DRef_wts'*(temp_uh_GQ.*temp_g_GQ)*GQ1DRef_wts;
                 term_A_elewise(temp_ele_id,1) = func_int;
                 term_A =  term_A +func_int;
                     
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
                 
                 term_B_elewise(temp_ele_id,1) = ACh_integral;
                 term_B = term_B+ACh_integral;
                     
             end
         end
     end
         
     
     if strcmp(func_type,'1') && strcmp(pde_ype,'1')
         Jh_star = term_A;
         ACh_star = term_B;
         Jh_star_elewise = term_A_elewise;
         ACh_star_elewise = term_B_elewise;
         
     elseif strcmp(func_type,'2') && strcmp(pde_ype,'1')
         Jh_star = term_B;
         ACh_star = term_A;
         Jh_star_elewise = term_B_elewise;
         ACh_star_elewise = term_A_elewise;
     else
         error('Wrong type')
     end
end
     
     
     
