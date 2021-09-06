function [grad_uH_x_GQ, grad_uH_y_GQ] ...
                    = GetGradUhstarFromUhstarGQ(uH_DG,GQ1DRef_pts,hx,hy)
                k = uH_DG.k;
                Lv = Vandermonde1D(k,GQ1DRef_pts);
                Dv = GradVandermonde1D(k,GQ1DRef_pts);
                NGQ = length(GQ1DRef_pts);
                [~,~,Nele] = size(uH_DG.dg_coeff);
                grad_uH_x_GQ = zeros(NGQ,NGQ,Nele);
                grad_uH_y_GQ = zeros(NGQ,NGQ,Nele);
                for ii = 1:uH_DG.Ny
                    for jj = 1: uH_DG.Nx
                        if uH_DG.mask(ii,jj) == 0
                            temp_ele_id = (ii-1)*uH_DG.Nx+jj;
                            temp_uh_coeff = uH_DG.dg_coeff(:,:,temp_ele_id);
                            temp_grad_uh_x = Dv * temp_uh_coeff * Lv' * 2.0/hx;
                            temp_grad_uh_y = Lv * temp_uh_coeff * Dv'* 2.0/hy;
                            grad_uH_x_GQ(:,:,temp_ele_id) = temp_grad_uh_x;
                            grad_uH_y_GQ(:,:,temp_ele_id) = temp_grad_uh_y;
                            
                        end
                    end
                end
                
end