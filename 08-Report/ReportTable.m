function ReportTable(varargin)
    fprintf('------------------------------\n')
    N = size(varargin{2});
    for ii = 1:2:nargin
        fprintf('%-9s  ',varargin{ii});
        if mod(ii,4) == 1
                fprintf('|');
        end
    end
    fprintf('\n');
    for nn = 1:N
        for ii = 2:2:nargin
            fprintf('%-9.2e  ',varargin{ii}(nn))
            if mod(ii,4) == 2
                fprintf('|');
            end
        end
        fprintf('\n')
    end
end