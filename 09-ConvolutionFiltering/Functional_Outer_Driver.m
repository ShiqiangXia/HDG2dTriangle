function [uh,qh,uhat,vh,ph,vhat,mymesh,Jh,ACh]=Functional_Outer_Driver(outer_mesh,...
        para, poly_order,...
        uH_star_inner,vH_star_inner,...
        Nadapt,err_inner,ndof_inner,area)
    pb_type = num2str(para.pb_type);
    [GQ1DRef_pts,GQ1DRef_wts] = GaussQuad(para.GQ_deg);
    
    posterior_estimate_method = 2;
    tol_adp = MyParaParse(para.extra_parameters,'tol_adp');
    % adaptive 
    mark_flag = 1; % 0: max marking strategy
                   % 1: bulk marking strategy Dorfler , 
                   % 2? equi distribution strategy
                   % 3: fraction marking strategy
    [uexact]=MyParaParse(para.pb_parameters,'uexact');
    
    mesh_list = zeros(Nadapt,1,numeric_t);
    
    err_list = zeros(Nadapt,1,numeric_t);
    
    tri_list  = zeros(Nadapt,1,numeric_t); % record # of triangles for each mesh
    cprintf('blue','adaptive steps \n');
    
    for ii = 1:Nadapt

        cprintf('blue','Mesh %d ... \n',ii)

        % Build initial mesh
        if ii == 1
            mymesh = outer_mesh;
        else
            if para.refine_flag <=0 % uniform refinement
                %% uniform refinement

                mymesh = mymesh.UniformRefine();

            else % refine based on marked elements
                %% adpvie refine
                r_f = GetRefineMethod(para.refine_flag); %'R', 'RG', 'NVB'
                mymesh = mymesh.Refine(marked_elements, r_f);

            end
        end


        flag_mesh_plot = 1;
        if flag_mesh_plot==1 && ii==Nadapt
            mymesh.Plot2(0,"Outer mesh " + num2str(ii));
%                 file_name = "k"+num2str(para.order)+"_outer_Mesh";
%                 if save_flag == 1
%                     savefig(gcf,file_name);
%                 end
        end


        mesh_list(ii) = GetDof(mymesh, poly_order);

        tri_list(ii) = mymesh.num_elements;
        if strcmp(pb_type(2),'0') % source problem
            %%  need to solve Primal and Adjoint two problems
            if strcmp(pb_type(3),'1') % Solve Poission source problem

                [source_f,uD,~]=MyParaParse(para.pb_parameters,'source_f','uD','uN');
                [source_g,vD,~]=MyParaParse(para.pb_parameters,'source_g','vD','vN');
                
                uD_inner = uH_star_inner;
                vD_inner = vH_star_inner;

                [uh,qh,uhat] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                poly_order, para.tau,source_f,uD,uD_inner); % temperorally set uD_inner to be uD

                [vh,ph,vhat] = HDG_SourcePbSolver_Elliptic(pb_type(3),mymesh,GQ1DRef_pts, GQ1DRef_wts,...
                poly_order, para.tau,source_g,vD,vD_inner);

            else
                error('Wrong problem type.')
            end

            [Jh,Jh_AC,ACh,ACh_elewise_list,Jh_elewise_list] = LinearFunctional_Elliptic(pb_type(4),pb_type(3),mymesh,...
                                                  uh,qh,uhat,source_f,...
                                                  vh,ph,vhat,source_g,...
                                                  GQ1DRef_pts,GQ1DRef_wts,poly_order,para.tau ,0);
            %{
            Jh_list(ii) = Jh;
            Jh_AC_list(ii) = Jh_AC;
            ACh_list(ii) = ACh;
            % Calculate error
            if err_cal_flag==1

                % error of uh, qh
                [uexact,qexact_1,qexact_2]=MyParaParse(para.pb_parameters,'uexact','qexact_1','qexact_2');

                [err_uh_list(ii),~] = L2Error_scalar(mymesh,uh,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    para.order,uexact);

                [err_qh_list(ii),~] = L2Error_vector(mymesh,qh,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    para.order,qexact_1,qexact_2);

                %PlotElementWiseValue(mymesh,err_uh_elewise,'err-uh elementwise' );

                [err_Jh_list(ii),err_Jh_AC_list(ii),err_Jh_elewise] ...
                    = Error_Functional(pb_type(4),para.pb_parameters,mymesh,GQ1DRef_pts,GQ1DRef_wts,Jh,Jh_AC,Jh_elewise_list);

                % error vh, ph;
                [vexact,pexact_1,pexact_2]=MyParaParse(para.pb_parameters,'vexact','pexact_1','pexact_2');
                [err_vh_list(ii),err_vh_elewise] = L2Error_scalar(mymesh,vh,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    para.order,vexact);
                %PlotElementWiseValue(mymesh,err_vh_elewise,'err-vh elementwise' );

            end
            %}


        end
        % HDG local post-processing for error estimation
        if para.post_process_flag == 1 

            if strcmp(pb_type(2),'0') && (strcmp(pb_type(4),'1') ||strcmp(pb_type(4),'2'))
                %%  poission source problem

                [uhstar,qhstar] = HDG_Local_Postprocess_Elliptic(mymesh,poly_order,uh,qh,uhat,para.tau,GQ1DRef_pts,GQ1DRef_wts);
                [vhstar,phstar] = HDG_Local_Postprocess_Elliptic(mymesh,poly_order,vh,ph,vhat,para.tau,GQ1DRef_pts,GQ1DRef_wts);
                %{ 
                err_uhstar_list(ii) = L2Error_scalar(mymesh,uhstar,...
                    GQ1DRef_pts,GQ1DRef_wts,1,...
                    para.order,uexact);

                err_qhstar_list(ii)= L2Error_vector(mymesh,qhstar,...
                   GQ1DRef_pts,GQ1DRef_wts,1,...
                    para.order,qexact_1,qexact_2);
                %}
                if posterior_estimate_method == 2
                    [est_terms_sum,est1,est2,est3,est4,est5,pos1,pos2,pos3,pos4]=...
                        Functional_Eh_Estimate_Residual_method(pb_type(4),pb_type(3),mymesh,...
                        uhstar,qhstar,source_f,...
                        vhstar,phstar,source_g,...
                        uh,qh,uhat,...
                        vh,ph,vhat,...
                        GQ1DRef_pts,GQ1DRef_wts,poly_order,para.tau,uexact);

                    %est_terms_sum_list(ii) = sum(est_terms_sum);
                end
            end
        end

        if para.post_process_flag == 1
            estimate_functinal_elewise = est_terms_sum+ACh_elewise_list;%+est_terms_sum ;%+%;
        else
            estimate_functinal_elewise = ACh_elewise_list;
        end

        % adaptive mesh refinement
        if para.refine_flag > 0
            percent = MyParaParse(para.extra_parameters,'percent');
            marked_elements = ACh_ErrEstimate(estimate_functinal_elewise,tol_adp,percent,mark_flag);
        end
        
        % evaluate the error
        [err_outer,~] = L2Error_scalar(mymesh,uh,...
                    GQ1DRef_pts,GQ1DRef_wts,0,...
                    poly_order,uexact);
        err_list(ii) = sqrt(err_outer^2 + err_inner^2);
        
        fprintf('err_inner: %.2e err_outer: %.2e\n',err_inner ,err_outer);
        
        
%         if mesh_list(ii)/(1-area)> 16*ndof_inner/area
%             mymesh.Plot2(0,"Outer mesh " + num2str(ii));
%             break
%         end        
        


    end
    
    figure;
    mesh_list = mesh_list + ndof_inner;
    plot(0.5*log10(mesh_list),log10(abs(err_list)),'--ro')
    legend('L2 error u-uh*')

    
end