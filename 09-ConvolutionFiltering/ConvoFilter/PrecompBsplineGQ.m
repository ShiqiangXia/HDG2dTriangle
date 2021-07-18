function AA = PrecompBsplineGQ(k,pts, GQ1DRef_pts)
    Npt = length(pts);
    NGQ = length(GQ1DRef_pts);
    J = 4*k+1;
    R = 2*k+1;
    AA = zeros(Npt,NGQ,R,J,numeric_t);
    
    for n = 1:Npt
        for rr = 1:R
            for jj = 1:J
                pts = 0.5*pts(n) - 0.5*GQ1DRef_pts - (jj-2*k-1) - (rr - k- 1);
                spline_pts =  Bspline(pts, k+1);
                AA(n,:,rr,J) = spline_pts';
            end
        end
    end
end