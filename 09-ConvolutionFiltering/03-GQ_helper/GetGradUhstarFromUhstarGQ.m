function [grad_uH_star_x_GQ, grad_uH_star_y_GQ] ...
                    = GetGradUhstarFromUhstarGQ(k,uH_star_GQ,GQ1DRef_pts,hx,hy)
                
                Lv = Vandermonde1D(k,GQ1DRef_pts);
                Iv = Lv\(V1D');
                IvT = Iv';
                Dv = GradVandermonde1D(k,GQ1DRef_pts);
end