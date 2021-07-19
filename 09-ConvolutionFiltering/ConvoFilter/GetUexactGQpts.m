function uexact_mat = GetUexactGQpts(uexact, GQ_x, GQ_y)
                
        [NGQ,num_ele] = size(GQ_x);
        uexact_mat = zeros(NGQ,NGQ,num_ele,numeric_t);
        for ii = 1:num_ele
            % get GQ points
            % total of NGQxNGQ points
            m = ones(NGQ,NGQ, numeric_t);
            x_list = m.*GQ_x(:,ii);
            x_list = reshape(x_list,[],1);

            y_list = m.*GQ_y(:,ii);
            y_list = reshape(y_list',[],1);

            uexact_pts = uexact([x_list,y_list]);

            uexact_pts = reshape(uexact_pts,[],NGQ);

            % make sure the matrix data matches with the structure of uh_GQ_pts
            % each row is the same y, x go from left to right
            uexact_pts = uexact_pts';
            uexact_mat(:,:,ii) = uexact_pts;
        end
end