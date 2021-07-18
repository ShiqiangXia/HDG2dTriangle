function out = Bspline(x,k)
% recursively compute Bspine of order k at x;
if (k == 1)
    out = zeros(size(x),numeric_t);
    idx = (x>=numeric_t('-0.5'))==(x<=numeric_t('0.5'));
    out(idx) = numeric_t('1');
else
    % compute recursively
    temp1 = Bspline(x+numeric_t('0.5'),k-1);
    temp2 = Bspline(x+numeric_t('-0.5'),k-1);
    s1 = k/numeric_t('2') + x;
    s2 = k/numeric_t('2') - x;
    out = numeric_t('1')/(k-1) * (s1.*temp1 + s2.*temp2);
end

end