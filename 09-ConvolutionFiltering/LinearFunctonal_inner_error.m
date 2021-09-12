function Eh_inner_elmtwise = LinearFunctonal_inner_error(...
                    qexact_1_GQ,qexact_2_GQ,qhstar_x_GQ,qhstar_y_GQ,...
                    pexact_1_GQ,pexact_2_GQ, phstar_x_GQ,phstar_y_GQ,...
                    Nx,Ny,mask,GQ1DRef_wts, hx, hy )
                
                scale_factor = hx * hy * numeric_t('1/4.0');
                Eh_inner_elmtwise = zeros(1, Nx*Ny,numeric_t);
                for ii = 1:Ny
                     for jj = 1: Nx
                         if mask(ii,jj) == 0
                             temp_ele_id = (ii-1)*Nx+jj;
                             formula = (qexact_1_GQ(:,:,temp_ele_id) - qhstar_x_GQ(:,:,temp_ele_id))...
                                 .*(pexact_1_GQ(:,:,temp_ele_id) - phstar_x_GQ(:,:,temp_ele_id)) ...
                                 + (qexact_2_GQ(:,:,temp_ele_id) - qhstar_y_GQ(:,:,temp_ele_id))...
                                 .*(pexact_2_GQ(:,:,temp_ele_id) - phstar_y_GQ(:,:,temp_ele_id));
                             Eh_inner_elmtwise(1,temp_ele_id) = scale_factor*GQ1DRef_wts'*(formula)*GQ1DRef_wts;
                             
                         end
                     end
                end

end
                