function RT_RHS = PostFlux_GetRHS(k,Jk,vertice_list,uhat_dir_list,A_qtau,Buuhat3,uh,qh,uhat,element_faces_list,tau)
    
    dimRT = (k+3)*(k+1);
    Nk_1 = (k*(k+1))/2;
    Nk = (k+1)*(k+2)/2;
    Nmu = k+1;
    Nuhat = k+1;
    
    dir_vec = GetDirVec(Nmu); % correct the uhat oritation 
    dir_vec = dir_vec';
    
    Id_mtrix = eye(Nmu,Nmu, numeric_t);
    
    RT_RHS = zeros(dimRT,1,numeric_t);
    
    RT_RHS(1:Nk_1) =Jk*A_qtau*qh(1:Nk);
    
    RT_RHS(Nk_1+1:2*Nk_1) = Jk*A_qtau*qh(Nk+1:end);
    
    [e_list,n1,n2,n3] = GetTriFaceInfo(vertice_list);
    for ii = 1:3
        
        start_id = 2*Nk_1 + (ii-1)*Nmu+1;
        end_id = start_id-1+Nmu;
        
        face_idx = element_faces_list(ii);
        temp_uhat = uhat((face_idx-1)*Nuhat+1:face_idx*Nuhat);
        
        if ii == 1
            norm_vec = n1;
        elseif ii == 2
            norm_vec = n2;
        else
            norm_vec = n3;
        end
        
        Bmu_q = 0.5*e_list(ii)* [norm_vec(1)*Buuhat3(:,:,ii);...
                norm_vec(2)*Buuhat3(:,:,ii)];
        Bmu_q = Bmu_q';
    
        Bmu_u = 0.5*e_list(ii)*tau*Buuhat3(:,:,ii);
        Bmu_u = Bmu_u';
        
        if uhat_dir_list(1,ii) == 0
            Bmu_uhat = 0.5*e_list(ii)*Id_mtrix.*dir_vec;
        else
            Bmu_uhat = 0.5*e_list(ii)*Id_mtrix;
        end
        
        RT_RHS(start_id:end_id) = Bmu_q*qh + Bmu_u * uh - tau*Bmu_uhat*temp_uhat ; 
        
        
    end
    
    
    
end