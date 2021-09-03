function [Jh_star, ACh_star] = LinearFunctional_uhstar_inner...
        ( uH_star_inner, source_f_GQ_pts,...
         vH_star_inner,source_g_GQ_pts,...
         GQ1DRef_wts, hx, hy)
     
     % compute two parts
     % part 1: Jh = (uh, g)
     % part 2: ACh = (f,vh)- (grad uh, grad vh)
     
     
     scale_factor = hx * hy * numeric_t('1/4.0');
     k = uH_star_inner.k;
     Lv = Vandermonde1D(k,GQ1DRef_wts);
     Dv = GradVandermonde1D(k,GQ1DRef_wts);
     Jh_star = 0;
     ACh_star = 0;
     
     for ii = 1:uH_star_inner.Ny
         for jj = 1: uH_star_inner.Nx
             if uH_star_inner.mask(ii,jj) == 0
                 temp_ele_id = (ii-1)*uH_star_inner.Nx+jj;
                 temp_uh_GQ = uH_star_inner.dg_gq_pts(:,:,temp_ele_id);
                 temp_g_GQ = source_g_GQ_pts(:,:,temp_ele_id);
                 
                 % (uh, g)
                 Jh_star =  Jh_star +...
                     scale_factor*GQ1DRef_wts'*(temp_uh_GQ.*temp_g_GQ)*GQ1DRef_wts;
                 
                 
                 temp_vh_GQ = vH_star_inner.dg_gq_pts(:,:,temp_ele_id);
                 temp_f_GQ = source_f_GQ_pts(:,:,temp_ele_id);
                 
                 temp_uh_coeff = uH_star_inner.dg_coeff(:,:,temp_ele_id);
                 temp_vh_coeff = vH_star_inner.dg_coeff(:,:,temp_ele_id);
                 temp_grad_uh_x = Dv * temp_uh_coeff * Lv';
                 temp_grad_uh_y = Lv * temp_uh_coeff * Dv';
                 temp_grad_vh_x = Dv * temp_vh_coeff * Lv';
                 temp_grad_vh_y = Lv * temp_vh_coeff * Dv';
                 
                 formula_grad_grad = temp_grad_uh_x.*temp_grad_vh_x ...
                     + temp_grad_uh_y.*temp_grad_vh_y ;
                 % (f,vh) - (grad uh, grad vh)
                 ACh_star = ACh_star  ...
                     + scale_factor*GQ1DRef_wts'*(temp_vh_GQ.*temp_f_GQ)*GQ1DRef_wts...
                     - scale_factor*GQ1DRef_wts'*(formula_grad_grad)*GQ1DRef_wts;
             end
         end
     end
         
     end
     
     
     
     
end