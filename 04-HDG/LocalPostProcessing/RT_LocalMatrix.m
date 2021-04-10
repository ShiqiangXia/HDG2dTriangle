function [A_qtau,RTuu,Buuhat3,RTuuhat3] = RT_LocalMatrix(k,GQ1DRef_pts,GQ1DRef_wts)
    
    % get the local matrix needed for the Raviart-Thomas space RT_k
    % A_qtau = (tau,q), where tau in P_{k-1}, q in RT_k
    % Buuhat3 = <u, mu>, u in P_k, mu in P_k(F)
    % RTuuhat3 = <u_ex, mu>, 
    % where RTuuhat3 has two components,
    % RTuuhat3(u_ex_id,mu_id,face_id,1) : u_ex = P_{i+1,k-i} 
    % RTuuhat3(u_ex_id,mu_id,face_id,1) : u_ex P_{i,k+1-i} for i= 0,1,..,k
    Nk = (k+1)*(k+2)/2;
    
    r_ind = Hierarchical_indx(k);

    temp_mat = eye(Nk,Nk,numeric_t);
    
    % (basis_P_k, basis_p_{k-1})
    A_qtau = temp_mat(r_ind,:);
    
    %RTuu = k+1 * ((k+1)*k/2) * 2
    % RTuu(:,:,1) = (r*P_{i,k-i}, basis_p_{k-1})
    %RTuu(:,:,2) =  (s*P_{i,k-i}, basis_p_{k-1})
    RTuu = Volume_Int_RTextrau_u(k,GQ1DRef_pts,GQ1DRef_wts);
    
    % <RT_basis, mu>
    
    Buuhat3 = Boundary_Int_u_uhat(k,GQ1DRef_pts,GQ1DRef_wts);
    
    
    RTuuhat3 = Boundary_Int_RTExtra_u_uhat(k,GQ1DRef_pts,GQ1DRef_wts);
    
    
    
    
    
    
    
    
    
end