function  ndof = GetInnerDof(mesh, outer_mesh,k )
     % interior dof
    INN = 11;
    Nu = (k+2)*(k+1)/2;
    Nq = 2*Nu;
    Nuhat = k+1;
    
    Nele = mesh.num_elements - outer_mesh.num_elements;
    Nface = mesh.num_faces - outer_mesh.num_faces + sum(outer_mesh.f_type==INN);
    ndof =  Nele*(Nu+Nq) + Nface*Nuhat;
    
%     % count elements
%     Nele_orderh = (Nx_coarse - 2*N_bd)*(Ny_coarse - 2*N_bd);
%     if N_corner_x >= N_bd && N_corner_y >= N_bd
%          Nele_order1 = (N_corner_x- N_bd)*(N_corner_y-N_bd);
%     else
%         Nele_order1 = 0;
%     end
%     
%     Nele = Nele_orderh - Nele_order1;
%     
%     % count faces
%     Nface_orderh = 2*Nele_orderh + (Nx_coarse - 2*N_bd) + (Ny_coarse - 2*N_bd);
%     Nface_order1 = 2*Nele_order1;
%     
%     Nface = Nface_orderh - Nface_order1 ;
%     
%     ndof =  Nele*(Nu+Nq) +Nface*Nuhat;
    
    
end