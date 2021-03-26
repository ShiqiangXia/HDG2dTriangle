function rlt = Project_F_to_Face(Jk,vertice_list,i,uhat_dir,edg_len,...
                    k,f,GQ1DRef_pts,GQ1DRef_wts)
                
     % project function f to the i-th face of the element defined by vertice_list
     rlt = zeros(k+1,1,numeric_t);
%      NGQ = length(GQ1DRef_pts);
     V1D = Vandermonde1D(k,GQ1DRef_pts);% legendre polynomail at GQ points
     
     [x_list,y_list] = GetFaceQaudPts(i, GQ1DRef_pts,Jk,vertice_list);
     f_VD = f([x_list,y_list]);
     
     for jj = 1:k+1
         rlt(jj,1) = 0.5*edg_len*GQ1DRef_wts'*(f_VD.*V1D(:,jj));
     end
     
     if uhat_dir == 0
         
         dir_vec = GetDirVec(k+1);
         
         rlt = rlt.*dir_vec;

     end
               
                
end