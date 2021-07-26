function out= Bspine_break_quadrature(deg_bspine,center_pt,x_l,x_r,r,w,Nu)
% compute integral using quadratures


% map to refer element
half_h = abs(x_r-x_l)*numeric_t('0.5');
mid    = (x_r+x_l) * numeric_t('0.5');
spline_pts = half_h * r + mid;
bspine_matrix = Bspline(spline_pts,deg_bspine);

% x_l <= center_pt - z/2 <= x_r
% z = (center_pt - spline_pts)*2

basis_pts = numeric_t('2')*(center_pt-spline_pts);
%basis_mtrix = my_vandermonde_u(basis_pts,Nu);
basis_mtrix = Vandermonde1D(Nu-1,basis_pts);


out = numeric_t('2')*half_h*bspine_matrix'*(w.*basis_mtrix);
end