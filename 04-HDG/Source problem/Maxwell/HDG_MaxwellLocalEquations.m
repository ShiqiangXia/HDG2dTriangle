function [Ns,M_Loc] = HDG_MaxwellLocalEquations(...
        Jk,vertice_list,tau_t,tau_n,mu,epsilon,omg...
        ,Aww,Awwr,Awws,Aww3,Bwuhat3,phat_dir_list)
    
    % Jk is the jacobian of the element
    % vertice_list [x1,y1;x2,y2;x3,y3]--?counterclockwise as the element vertice index
    % HDG tau
    
    Nw = size(Aww,1); 
    Np= Nw;
    Nu = 2*Nw;
    N_local = Nw+Nu+Np;
    
    Nface = size(Bwuhat3,2);
    N_global = 2*Nface;
     
    dir_vec = GetDirVec(Nface); % correct the uhat oritation 
    dir_vec = dir_vec';
    
    % ph_hat has to have a unqiue oritation, we use the oritentiaon:
    % small_verice_idx ---> large_vertice_idx
    % the local matrix Buuhat3 = <phat,p> is using the counterclock wise
    % orientation. it's not unique in two triangles, so we need to fix it. 
    % the difference is legendre_i(-x) = (-1)^i*legendre_i(x)
    
    %---------------------------------------------------------------------
    % M_loc
    % (w_j,w_i)
    Mww = Jk*mu*Aww;
    
    %(u_j,curl x w_i)
    % get the affine map and it's inverse
    
    [~,Inv_AffineMap] = Ref_Tri_Map(Jk,vertice_list);
    r_x = Inv_AffineMap(1,1);
    r_y = Inv_AffineMap(1,2);
    s_x = Inv_AffineMap(2,1);
    s_y = Inv_AffineMap(2,2);
                        
    temp1 =  Jk * (Awwr*r_y + Awws*s_y);
    temp2 = -Jk * (Awwr*r_x + Awws*s_x);
    
    Muw = [temp1;temp2];
    
    
    % <tau_t u_jxn, u_ixn>
    [e_list,n1,n2,n3] = GetTriFaceInfo(vertice_list);
    
    Mu1u1 = 0.5*tau_t*(...
          e_list(1)*Aww3(:,:,1)*n1(2)*n1(2)...
        + e_list(2)*Aww3(:,:,2)*n2(2)*n2(2)...
        + e_list(3)*Aww3(:,:,3)*n3(2)*n3(2));
    
    Mu1u2 = - 0.5*tau_t*(...
          e_list(1)*Aww3(:,:,1)*n1(2)*n1(1)...
        + e_list(2)*Aww3(:,:,2)*n2(2)*n2(1)...
        + e_list(3)*Aww3(:,:,3)*n3(2)*n3(1));
%     Mu2u1 = - 0.5*tau_t(...
%           e_list(1)*Aww3(:,:,1)*n1(1)*n1(2)...
%         + e_list(2)*Aww3(:,:,2)*n2(1)*n2(2)...
%         + e_list(3)*Aww3(:,:,3)*n3(1)*n3(2));

    Mu2u1 = Mu1u2;
    
    Mu2u2 = 0.5*tau_t*(...
          e_list(1)*Aww3(:,:,1)*n1(1)*n1(1)...
        + e_list(2)*Aww3(:,:,2)*n2(1)*n2(1)...
        + e_list(3)*Aww3(:,:,3)*n3(1)*n3(1));
    
    Muu = [Mu1u1,Mu1u2;Mu2u1,Mu2u2] + epsilon * omg^2 * Jk* blkdiag(Aww,Aww);
    
    temp_3 = epsilon*Jk*(Awwr*r_x + Awws*s_x);
    temp_4 = epsilon*Jk*(Awwr*r_y + Awws*s_y);
    Mpu = [temp_3,temp_4];
    
    Mpp = tau_n*epsilon*0.5*(...
        e_list(1)*Aww3(:,:,1)...
        + e_list(2)*Aww3(:,:,2)...
        + e_list(3)*Aww3(:,:,3));
    
    zero_mat = zeros(Nw,Nw,numeric_t);
    
    M_Loc =[Mww, -Muw', zero_mat;...
            Muw,  Muu,    -Mpu' ;...
        zero_mat, Mpu,     Mpp ];
        
    
    
    
    
    %---------------------------------------------------------------------
    % N_1,N_2,N_3
    Ns = zeros(N_local,N_global,3,numeric_t);
    for tt = 1:3
        
        if tt ==1
            n_vec = n1;
        elseif tt == 2
            n_vec = n2;
        else
            n_vec = n3;
        end
        
        Z_wuhat_t = Bwuhat3(:,:,tt);
        
        Z_uuhat_t = tau_t*[Bwuhat3(:,:,tt)*n_vec(2);-Bwuhat3(:,:,tt)*n_vec(1)];
        
        Z_uphat = -epsilon*[Bwuhat3(:,:,tt)*n_vec(1);Bwuhat3(:,:,tt)*n_vec(2)];
        
        Z_pphat = tau_n*epsilon*Bwuhat3(:,:,tt);
        if phat_dir_list(1,tt) == 0  % reverse uhat direction
            Z_wuhat_t = Z_wuhat_t.*dir_vec;
            Z_uuhat_t = Z_uuhat_t.*dir_vec;
            Z_uphat = Z_uphat.*dir_vec;
            Z_pphat = Z_pphat.*dir_vec;
            
            Z_wuhat_t = - Z_wuhat_t; % delta_Th
            Z_uuhat_t = - Z_uuhat_t;
        end
        zero_mat2 = zeros(Nw,Nface,numeric_t);
        
        N = [Z_wuhat_t,  zero_mat2;...
             Z_uuhat_t,  Z_uphat;...
             zero_mat2,  Z_pphat];
         
        
        Ns(:,:,tt) = 0.5*e_list(tt)*N;
        
    end
    

    
end