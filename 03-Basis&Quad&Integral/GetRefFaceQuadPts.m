function [r_list,s_list] = GetRefFaceQuadPts(face_local_id,GQ1DRef_pts)
    
    NGQ = length(GQ1DRef_pts);
    if face_local_id == 1
     % edge 1:  r=[-1,1], s= -1
        r_list = GQ1DRef_pts;
        s_list = -1*ones(NGQ,1,numeric_t);

     elseif face_local_id == 2
     % edge 2:  r=-2, s= [-1,1]
     
        s_list = GQ1DRef_pts;
        r_list = -s_list;
        
   
     elseif face_local_id == 3
     % edge 3:  r=-1, s= [1,-1]  
        r_list = -1*ones(NGQ,1,numeric_t);
        s_list = -GQ1DRef_pts;
  
     else 
         error('Wrong face index');
    end
    
end
