function rlt = Project_F_to_Face(Jk,vertice_list,i,...
                    k,f,GQ1DRef_pts,GQ1DRef_wts)
                
     % project function f to the i-th face of the element defined by vertice_list
     rlt = zeros(k+1,1,numeric_t);
     NGQ = length(GQ1DRef_pts);
     V1D = Vandermonde1D(k,GQ1DRef_pts);% legendre polynomail at GQ points
     [e_list,~] = GetTriFaceInfo(vertice_list);
     edg_len = e_list(i);
     
     if i == 1
     % edge 1:  r=[-1,1], s= -1
        r_list = GQ1DRef_pts;
        s_list = -1*ones(NGQ,1,numeric_t);
        [x_list,y_list] = RStoXY(r_list,s_list,Jk,vertice_list);
        f_VD = f([x_list,y_list]);
        
     elseif i == 2
     % edge 2:  r=[-1,1], s= -r
        r_list = GQ1DRef_pts;
        s_list = -r_list;
        [x_list,y_list] = RStoXY(r_list,s_list,Jk,vertice_list);
        f_VD = f([x_list,y_list]);
     elseif i == 3
     % edge 3:  r=-1, s= [-1,1]  
        r_list = -1*ones(NGQ,1,numeric_t);
        s_list = GQ1DRef_pts;
        [x_list,y_list] = RStoXY(r_list,s_list,Jk,vertice_list);
        f_VD = f([x_list,y_list]);
     else 
         error('Wrong face index');
     end
     
     for jj = 1:k+1
         rlt(jj,1) = 0.5*edg_len*GQ1DRef_wts'*(f_VD.*V1D(:,jj));
     end
               
                
end