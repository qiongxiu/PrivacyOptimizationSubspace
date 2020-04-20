%%%plot the comparison of admm and pdmm under different privacy levels 
function subPert_privplot(error_th,ave_Res,name)
intVal1=ave_Res.intVal(1);
intVal2=ave_Res.intVal(2);
intVal3=ave_Res.intVal(3);
N=10;
Markers = {'o','x','s','v','d','^','>','<'};
figure; 
set(gca,'fontsize',15)
for i=1:3
    plot( ave_Res.admm(i).transmission(:),ave_Res.admm(i).MSE_error(:),strcat('b-'),'Marker',Markers{i},'MarkerIndices',1:intVal1:length(ave_Res.admm(i).MSE_error),'MarkerSize',10,'linewidth',1.1);
        hold on 
end 
ylim([error_th inf ])
% xlim([0 60])
set(gca,'fontsize',15)
grid on 
legend({'p-ADMM: $7*{10}^{-7}$','p-ADMM: $7*{10}^{-5}$','p-ADMM: $7*{10}^{-3}$'},'location','northeast','FontSize',15,'Interpreter','Latex')
set(gca, 'YScale', 'log')
xlabel ('Transmissions'); ylabel ('||x^{(k)}-x*||^{2}')
title(strcat('Privacy Levels: Distributed ', name))
set(gca, 'FontSize', 12)
set(gca,'linewidth',2);

figure; 
set(gca,'fontsize',15)
for i=1:3

    plot( ave_Res.pdmm(i).transmission(:),ave_Res.pdmm(i).MSE_error(:),strcat('b-'),'Marker',Markers{i},'MarkerIndices',1:intVal2:length(ave_Res.pdmm(i).MSE_error),'MarkerSize',10,'linewidth',1.1);
        hold on 
end
set(gca,'fontsize',15)
ylim([error_th inf])
xlim([sum(sum(ave_Res.Geograph.weight)) inf])
grid on 
legend({'p-PDMM: $7*{10}^{-7}$','p-PDMM: $7*{10}^{-5}$','p-PDMM: $7*{10}^{-3}$'},'location','northeast','FontSize',15,'Interpreter','Latex')
set(gca, 'YScale', 'log')
xlabel ('Transmissions'); ylabel ('||x^{(k)}-x*||^{2}')
title(strcat('Privacy Levels: Distributed ', name))
set(gca, 'FontSize', 12)
set(gca,'linewidth',2);

figure; 
set(gca,'fontsize',15)
for i=1:3  
    plot( ave_Res.dual(i).transmission(:),ave_Res.dual(i).MSE_error(:),strcat('b-'),'Marker',Markers{i},'MarkerIndices',1:intVal3:length(ave_Res.dual(i).MSE_error),'MarkerSize',10,'linewidth',1.1);
    hold on 
end 
set(gca,'fontsize',15)
ylim([error_th inf])
% xlim([0 60])
grid on 
legend({'p-Dual: $7*{10}^{-7}$','p-Dual: $7*{10}^{-5}$','p-Dual: $7*{10}^{-3}$'},'location','northeast','FontSize',15,'Interpreter','Latex')
set(gca, 'YScale', 'log')
xlabel ('Transmissions'); ylabel ('||x^{(k)}-x*||_{2}^{2}')
title(strcat('Privacy Levels: Distributed ', name))
set(gca, 'FontSize', 12)
set(gca,'linewidth',2);

