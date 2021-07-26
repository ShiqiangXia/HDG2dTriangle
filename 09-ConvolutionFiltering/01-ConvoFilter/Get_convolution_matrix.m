function Conv_matrix=Get_convolution_matrix(pts,kk,Nu,GQ_pts,GQ_weights)


AA = Convolution_matrix(kk,Nu,pts,GQ_pts,GQ_weights);


Cr1 = numeric_t('[-(1/12), 7/6, -(1/12)]');
Cr2 = numeric_t('[37/1920, -(97/480), 437/320, -(97/480), 37/1920]');
Cr3 = numeric_t('[-(41/7560), 311/5040, -(919/2520), 12223/7560, -(919/2520), 311/5040,-(41/7560)]');
Cr4 = numeric_t('[153617/92897280, -(35411/1658880), 3153959/23224320, -(6803459/11612160), 18017975/9289728, -(6803459/11612160),3153959/23224320, -(35411/1658880), 153617/92897280]');
Cr5 = numeric_t('[-(4201/7983360), 30773/3991680, -(20813/380160), 2825/11088,-(1179649/1330560), 1569217/665280, -(1179649/1330560), 2825/11088,-(20813/380160), 30773/3991680, -(4201/7983360)]');

if kk == 1
    cr = Cr1;
elseif kk ==2
    cr = Cr2;
elseif kk == 3
    cr = Cr3;
elseif kk == 4
    cr = Cr4;
elseif kk == 5
    cr = Cr5;
else
    error("Order %d B-slpine kernel is not implemented yet.\n",kk);
end
%str = ['c', num2str(kk)];
%cr = numeric_t(struct2array( load("Convolution_Cr_deg1_5.mat",str)));

Conv_matrix = zeros(length(pts),Nu,4*kk+1,numeric_t);

for ss = 1:2*kk+1
    Conv_matrix(:,:,:) = Conv_matrix(:,:,:) + cr(ss)*squeeze(AA(:,:,ss,:))*numeric_t('0.5');
end

end