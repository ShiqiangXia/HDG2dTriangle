function [rlt,Cs] = Post_Estimate(dofs, Jhs, Dhs)
% Idea: assume J-Jh = CDh, use two data points to estimate C and then get a better approximation C*Dh 

   % step 1: get the type for data 2: end-1 for Err_Jh and Dh
    %  type 1() \\      type 2 (): \/   type 3: /\    type 4: // (not seen)
    %                          
    % step 2: Compare types
    % if type match, assume paranell
    % if not:
    %       deal with each case
    
    method_flag = 2; % 2
    N = size(dofs,1);
    rlt = zeros(N,1,numeric_t);
    Cs = zeros(N,1,numeric_t);
    if N == 1
        return
    end
    
    dff_Jh = Jhs(1:end-1) - Jhs(2:end);
    if method_flag == 1
        %% 
        for ii = 1:N
            if ii>1 && ii<N-1

                Jh_type = CheckType(ii,abs(dff_Jh));
                Dh_type = CheckType(ii,abs(Dhs));

                if Jh_type == 1 && Dh_type ==2

                    x1 = 0.5*log10(dofs(ii-1));
                    y1 = log10(abs(Dhs(ii-1)));

                    x3 = 0.5*log10(dofs(ii+1));
                    y3 = log10(abs(Dhs(ii+1)));

                    slope = (y1 - y3)/(x1-x3);

                    x2 = 0.5*log10(dofs(ii));
                    y2 = slope*(x2-x1)+y1;

                    Cs(ii,1) = (Jhs(ii) - Jhs(ii-1))/(10^y1 - 10^y2);

                    rlt(ii,1) = Cs(ii,1)*(10^y2);

                elseif Jh_type == 2 && Dh_type ==1

                     Cs(ii,1) = (Jhs(ii) - Jhs(ii-1))/(Dhs(ii-1) - Dhs(ii));
                     rlt(ii,1) = Cs(ii,1)*Dhs(ii); 

                elseif Jh_type == 2 && Dh_type ==3
                    x1 = 0.5*log10(dofs(ii-1));
                    y1 = log10(abs(Dhs(ii-1)));

                    x3 = 0.5*log10(dofs(ii+1));
                    y3 = log10(abs(Dhs(ii+1)));

                    slope = (y1 - y3)/(x1-x3);


                    x2 = 0.5*log10(dofs(ii));
                    y2 = slope*(x2-x1)+y1;

                    Cs(ii,1) = (Jhs(ii-1) - Jhs(ii+1))/(Dhs(ii+1) - Dhs(ii-1));

                    rlt(ii,1) = Cs(ii,1)*(10^y2);
                else

                    Cs(ii,1) = (Jhs(ii) - Jhs(ii+1))/(Dhs(ii+1) - Dhs(ii));
                    rlt(ii,1) = Cs(ii,1)*Dhs(ii);

                end

            elseif ii == 1 || ii == N-1
                Cs(ii,1) = (Jhs(ii) - Jhs(ii+1))/(Dhs(ii+1) - Dhs(ii));
                rlt(ii,1) = Cs(ii,1)*Dhs(ii);
            else
                Cs(ii,1) = (Jhs(ii) - Jhs(ii-1))/(Dhs(ii-1) - Dhs(ii));
                rlt(ii,1) = Cs(ii,1)*Dhs(ii);

            end
        end
    elseif method_flag == 2
        %% 
        for ii = 1:N
            if ii>1 && ii<N

                Dh_type = CheckType(ii,abs(Dhs));

                if Dh_type ==2

                    x1 = 0.5*log10(dofs(ii-1));
                    y1 = log10(abs(Dhs(ii-1)));

                    x3 = 0.5*log10(dofs(ii+1));
                    y3 = log10(abs(Dhs(ii+1)));

                    slope = (y1 - y3)/(x1-x3);

                    x2 = 0.5*log10(dofs(ii));
                    y2 = slope*(x2-x1)+y1;

                    Cs(ii,1) = (Jhs(ii) - Jhs(ii-1))/(10^y1 - 10^y2);

                    rlt(ii,1) = Cs(ii,1)*(10^y2);


                elseif Dh_type ==3
                    x1 = 0.5*log10(dofs(ii-1));
                    y1 = log10(abs(Dhs(ii-1)));

                    x3 = 0.5*log10(dofs(ii+1));
                    y3 = log10(abs(Dhs(ii+1)));

                    slope = (y1 - y3)/(x1-x3);


                    x2 = 0.5*log10(dofs(ii));
                    y2 = slope*(x2-x1)+y1;

                   
                    Cs(ii,1) = (Jhs(ii-1) - Jhs(ii+1))/(Dhs(ii+1) - Dhs(ii-1));

                    rlt(ii,1) = Cs(ii,1)*(10^y2);
                else

                    Cs(ii,1) = (Jhs(ii) - Jhs(ii+1))/(Dhs(ii+1) - Dhs(ii));
                    rlt(ii,1) = Cs(ii,1)*Dhs(ii);

                end

            elseif ii == 1 
                Cs(ii,1) = (Jhs(ii) - Jhs(ii+1))/(Dhs(ii+1) - Dhs(ii));
                rlt(ii,1) = Cs(ii,1)*Dhs(ii);
            else
                Cs(ii,1) = (Jhs(ii) - Jhs(ii-1))/(Dhs(ii-1) - Dhs(ii));
                rlt(ii,1) = Cs(ii,1)*Dhs(ii);

            end
        end

    end
 
end