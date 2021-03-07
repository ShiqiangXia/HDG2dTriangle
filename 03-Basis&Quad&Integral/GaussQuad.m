function [x,w] = GaussQuad(N)
    
% Get the Gauss Quadrature with N points and weights
% accuray  2N-1

    [x,w] = JacobiGQ(0,0,N-1);

end