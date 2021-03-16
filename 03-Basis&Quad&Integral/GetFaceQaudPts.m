function [x_list,y_list] = GetFaceQaudPts(face_local_id, GQ1DRef_pts,Jk,vertice_list)
    
    [r_list,s_list] = GetRefFaceQuadPts(face_local_id,GQ1DRef_pts);
    
    [x_list,y_list] = RStoXY(r_list,s_list,Jk,vertice_list);
     
end