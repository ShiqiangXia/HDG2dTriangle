function VisulaizeElewiseError(...
                    data0,data1,data2,...
                    mask,Nx,Ny)
                
       for ii = 1:Ny
         for jj = 1: Nx
             if mask(ii,jj) == 0
                 ele_id = (ii-1)*Nx+jj;
                 data1(ele_id,1) = data1(ele_id,1) ...
                     - (data0(ele_id,1) + data0(ele_id + Nx*Ny,1));
                
             end
         end
       end
       
       figure;
       plot(1:Nx*Ny,data1,'--rx');
       legend('J-Jh_star_inner')
end