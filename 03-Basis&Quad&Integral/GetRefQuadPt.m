function [a_list,b_list,Jacobian_rs_to_ab]= GetRefQuadPt(GQ1DRef_pts)
    NGQ = length(GQ1DRef_pts);
    m = ones(NGQ,NGQ, numeric_t);
    m = m.*GQ1DRef_pts;
    a_list = reshape(m',[],1);
    b_list = reshape(m,[],1);
    
    Jacobian_rs_to_ab = (1-b_list)/2;
    Jacobian_rs_to_ab = reshape(Jacobian_rs_to_ab,[],NGQ);
end