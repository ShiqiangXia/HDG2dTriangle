function uh_coeff = GetUhL2ProjectionCoarseMesh(poly_k, coarse, fine,...
        uh, GQ1DRef_pts, GQ1DRef_wts, hx, hy)
    % Goal: L2 projection of uh defined on fine mesh (triangle) to coarse mesh(square)
    % Output: uh_coeff is a matrix of Nk x Nk x num_element
    % (i,j) is for the basis P_i(x)*P_j(y)
    % so each row is the same x basis, diff y basis
    % WARNING: right now the code ONLY works for unit square domain
    % ASSUME THE COARSE AND FINE MESH IS THE SAME MESH
    
    % Get relation mat
    %relation_mat = BuildRelationof2Meshes(coarse, fine);
    
    
    num_coarse = coarse.num_elements ;
    % this follwong only work for square domain
    Nsquare = num_coarse / 2; % two triangels make a square
    Ns = sqrt(Nsquare);
    Nk = poly_k+1;
    NGQ = length(GQ1DRef_pts);
    scale_factor = 4/(hx*hy);
    
    [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts);
    V2D = Vandermonde2D(poly_k,a_list,b_list);  % scalar Pk
    [r_list,s_list] = ABtoRS(a_list,b_list);
    
    uh_coeff = zeros(Nk,Nk,Nsquare, numeric_t);
    
    for ee = 1:Ns
        for rr = 1:Ns
            square_idx = (ee-1)*Ns + rr;
            xmid = (rr-1)*hx + 0.5*hx;
            ymid = (ee-1)*hy + 0.5*hy;
            % do (uh_proj, basis_ij) = sum (uh, basis_ij)_K for all elements
            % K contained in this square
%             tri_elements = [relation_mat{square_idx},...
%                 relation_mat{square_idx + Nsquare}];

            tri_elements = [square_idx,...
                square_idx + Nsquare];
            
            n = length(tri_elements);
            
            for tt = 1:n
                % get Gauss points on this triangle
                element_idx = tri_elements(tt);
                uh_pts = V2D * (uh(:,element_idx));
                uh_pts = reshape(uh_pts,[],NGQ);
                
                temp_element = fine.element_list(element_idx,:);
                vertice_list = fine.vertices_list(temp_element(:),:);
                Jk = fine.Jacobian_list(element_idx);
                [x_phy_list,y_phy_list] = RStoXY(r_list,s_list,Jk,vertice_list);
                %labs = 1:NGQ*NGQ;
                %plot3(x_phy_list,y_phy_list, 0 *labs,'r*')
                
                x_list = (x_phy_list - xmid)/(hx*0.5); % map the physical GQ points on triangle to [-1,1]
                y_list = (y_phy_list - ymid)/(hy*0.5);
                
                V1D_x = Vandermonde1D(poly_k,x_list);
                V1D_y = Vandermonde1D(poly_k,y_list);
                
                mat = zeros(Nk,Nk, numeric_t);
                for ii = 1:Nk
                    for jj = 1:Nk
                        basis_ij = V1D_x(:,ii) .* V1D_y(:,jj);
                        basis_ij = reshape(basis_ij,[],NGQ);

                        temp = (uh_pts.*basis_ij);
                        mat(ii,jj) = Jk*GQ1DRef_wts'*(temp.*Jacobian_rs_to_ab)*GQ1DRef_wts;

                    end
                end

                uh_coeff(:,:,square_idx) = uh_coeff(:,:,square_idx) + mat;
            end
            
        end
    end
    
    uh_coeff = uh_coeff * scale_factor;
    
end