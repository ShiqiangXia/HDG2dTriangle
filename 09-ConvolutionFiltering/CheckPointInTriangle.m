function t = CheckPointInTriangle(P, P1,P2,P3)
    % ref
    % % https://www.mathworks.com/matlabcentral/answers/277984-check-points-i%C3%A5nside-triangle-or-on-edge-with-example
    s = det([P1-P2;P3-P1]);
    t = s*det([P3-P;P2-P3])>=0 ...
        & s*det([P1-P;P3-P1])>=0 ...
        & s*det([P2-P;P1-P2])>=0;
end