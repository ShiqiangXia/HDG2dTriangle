% Oliver Xia 01/08/2020
% Return Convolution matrix 
% #points x Nu x (2k+1) x (4k+1)

function AA = Convolution_matrix(k,Nu,pts,GQ_points,GQ_weights)

% use Bspine of order k+1;
% Polynomial degree is Nu-1;
% pts are the points of interest;
% pts need to be a COLUMN vector

%--------------------------------------------------------------------------
%             Implementation details
%
% Notice that Bspine_k+1 is a piecewise polynomial with nodes
%
% I = [-(k+1)/2, -(k-1)/2, ...., (k-1)/2, (k+1)/2]
%   = [x_0, x_1, ......x_{k+1}]
%  k+1 intervals with width 1

% Bspine_k+1 is only in C_{k-1}
% When we compute Int(Bspline_k+1,basis_u), we need to use the Gauss
% Quadratre for each interval to get the right accuracy
%
% Therefor for given r,j,point pt_n
% center_point s:= pt_n/2 - j -r
% intergral domain [s-1/2,s+1/2]
% this domain is located somewhere in the intervales
%
% case 1 [s-1/2,s+1/2] is outside --> integral = 0
% case 2 [s-1/2,s+1/2] is part inside --->integral [x_0,s+1/2] or [s-1/2,x_{k+1}]
% case 3 [s-1/2,s+1/2] is inside but not overlap with nodes
%        split integral : [s-1/2, x_i] and [x_i,s+1/2]
% case 3 [s-1/2,s+1/2]  happens to be nodes --> integral [x_i,x_{i+1}]
%
% more details see notes
%--------------------------------------------------------------------------
deg_bspine = k+1;


n = length(pts);
AA = zeros(n,Nu,2*k+1,4*k+1,numeric_t);



for rr = -k:1:k
    for jj = -2*k:1:2*k
        for nn = 1:n
        % compute pt_n/2 - jj -r
            center_pt = pts(nn)/numeric_t('2') - jj - rr;
            %[breaks,spline_flag] = Bspline_breaks(deg_bspine,temp);
        
        % compute Bspine_k+1 (pt_n/2 - jj - r_n/2 -r)
        % matrix: each row same point, diffent Gauss points rn
         
             
            bspine_matrix = Bspline_int(deg_bspine,center_pt,Nu,GQ_points,GQ_weights);
            %bspine_matrix = Bspine(temp,deg_bspine);
            
            
        
        % numerical quadratue
            AA(nn,:,rr+k+1,jj+2*k+1)=bspine_matrix;
        end
    end
end

end