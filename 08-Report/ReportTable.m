function ReportTable(varargin)
    fprintf('------------------------------\n')
    N = size(varargin{2});
    for ii = 1:2:nargin
        fprintf('%-10s  ',varargin{ii});
    end
    fprintf('\n');
    for nn = 1:N
        for ii = 2:2:nargin
            fprintf('%-10.2e  ',varargin{ii}(nn))
        end
        fprintf('\n')
    end
end