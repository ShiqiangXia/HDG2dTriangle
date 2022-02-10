

hold on
% k = 1
% mesh_list = [3.68e+02,1.50e+03,6.08e+03,2.44e+04];
% err_lamh_AC_Dh_list = [3.99e-03,1.66e-03,1.18e-04,7.38e-06];

% % k=2
% mesh_list =[6.96e+02,2.83e+03,1.14e+04,4.59e+04  ];
% err_lamh_AC_Dh_list =[6.21e-04,8.91e-06,1.35e-07 ,2.08e-09];

% k= 3
mesh_list = [1.12e+03,4.54e+03,1.83e+04,7.35e+04];
err_lamh_AC_Dh_list = [ 4.94e-06,1.86e-08,7.12e-11,3.43e-13];



plot(0.5*log10(mesh_list),log10(abs(err_lamh_AC_Dh_list)),'-m*',...
                            'MarkerSize',10,'LineWidth',1);
                        
legend('$|\lambda-\lambda_h|$',...
        '$|\mathcal{E}_h|$',...
        '$\mathcal{E}_{h,|\cdot|}$',...
        '$|\lambda-\lambda_h-\mathcal{E}_h|$',...
        'Interpreter','latex','FontSize',15)                       