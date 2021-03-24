function eig_err=ParseEigenError(mesh_list,err_lamh_list,tag_text,varargin)
    [N,Neig] = size(err_lamh_list);
    
    eig_err = cell(1,4*Neig);
    text_order = 'order';
    for ii = 1:Neig
        temp_tag = append(tag_text,num2str(ii));
        temp_eig_err = err_lamh_list(:,ii);
        if nargin<=3
            temp_order = GetOrder(mesh_list,temp_eig_err);
        else
            temp_order = zeros(N,1);
        end
        eig_err{1,(ii-1)*4+1} = temp_tag;
        eig_err{1,(ii-1)*4+2} = temp_eig_err;
        eig_err{1,(ii-1)*4+3} = text_order;
        eig_err{1,(ii-1)*4+4} = temp_order;
        
        
    end
    
end