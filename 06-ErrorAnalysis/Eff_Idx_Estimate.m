function Eff_Idx_Estimate(dofs, Jhs, Dhs)
    figure;
    plot(log10(dofs),log10(abs(Dhs)),'--rx',...
        log10(dofs(1:end-1)),log10(abs(Jhs(1:end-1)-Jhs(2:end))),'--bo')
    
    % step 1: get the type for data 2: end-1 for Err_Jh and Dh
    %  type 1() \\      type 2 (): \/   type 3: /\    type 4: // (not seen)
    %                          
    % step 2: Compare types
    % if type match, assume paranell
    % if not:
    %       deal with each case
end