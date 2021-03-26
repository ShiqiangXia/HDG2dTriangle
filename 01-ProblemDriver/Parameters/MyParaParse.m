function  varargout= MyParaParse(para,varargin)
    % Parse pameters from para based on index in varagin
    % para: parmeters to parse, format: pair 'index', value
    % varargin: index
    %
    % Example: para = {'x1',5,'zz',70,'x3',6,'y',[1,2], 'x2','somestring'}
    % MyParaParse(para,'x1','x2','x3','y') will return {5,'somestring',6,[1,2]}
    if nargin <2
        error('Inputs of MyParaPrase is incorrect.At least two inputs are needed. \n')
    else
        p = inputParser;
        for j = 1:nargin-1
            addParameter(p,varargin{j},'error');
        end
        p.KeepUnmatched = true;
        parse(p,para{:})
        
       
        for j = 1:nargin-1
            rlt = p.Results.(varargin{j});
            if strcmp(rlt,'error')
                error('Index %s is not find in para',varargin{j})
            else
                varargout{j} = rlt;
            end
        end
        
    end
     
end