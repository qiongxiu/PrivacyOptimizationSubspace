%%%plot the convergence of admm, pdmm and dual ascent under different
%%%applications
function subPert_converplot(error_th,aveRes,Pmax,name)
output_a1=aveRes.output_a1;
output_a2=aveRes.output_a2;
output_p1=aveRes.output_p1;
output_p2=aveRes.output_p2;
output_d1=aveRes.output_d1;
output_d2=aveRes.output_d2;
intVal1=aveRes.intVal(1);
intVal2=aveRes.intVal(2);
intVal3=aveRes.intVal(3);

N=10;
Markers = {'+','o','*','x','v','d','^','s','>','<'};
figure; 
set(gca,'fontsize',15)
p= plot(output_a1.transmission(:),output_a1.MSE_error(:),strcat('b-'),'Marker',Markers{3},'MarkerIndices',1:intVal1:length(output_a1.MSE_error),'MarkerSize',10,'linewidth',1.1);
hold on; plot(output_a1.transmission(:),output_a1.Z_Con_error(:),strcat('b-'),'Marker',Markers{6},'MarkerIndices',1:intVal1:length(output_a1.Z_Con_error(:)),'MarkerSize',10,'linewidth',1.1);
hold on; plot(output_a1.transmission(:),output_a1.Z_nCon_error(:),strcat('b-'),'Marker',Markers{9},'MarkerIndices',1:intVal1:length(output_a1.Z_nCon_error(:)),'MarkerSize',10,'linewidth',1.1);
hold on; plot(output_a2.transmission(:),output_a2.MSE_error(:),strcat('r-'),'Marker',Markers{3},'MarkerIndices',1:intVal1:length(output_a2.MSE_error(:)),'MarkerSize',10,'linewidth',1.1); 
hold on; plot(output_a2.transmission(:),output_a2.Z_Con_error(:),strcat('r-'),'Marker',Markers{6},'MarkerIndices',1:intVal1:length(output_a2.Z_Con_error(:)),'MarkerSize',10,'linewidth',1.1);
hold on; plot(output_a2.transmission(:),output_a2.Z_nCon_error(:),strcat('r-'),'Marker',Markers{9},'MarkerIndices',1:intVal1:length(output_a2.Z_nCon_error(:)),'MarkerSize',10,'linewidth',1.1);
ylim([error_th Pmax^2*10])
% xlim([0 60])
grid on 
legend({'Proposed: x^{(k)}-x*','Proposed: (\Pi_{H})v^{(k)}-v*','Proposed: (I-\Pi_{H})v^{(k)}','Non-private: x^{(k)}-x*','Non-private: (\Pi_{H})v^{(k)}-v*','Non-private: (I-\Pi_{H})v^{(k)}'},'location','northeast','FontSize',10)
set(gca, 'YScale', 'log')
xlabel ('Transmissions'); ylabel ('MSE')
title(strcat('Distributed ', name,' of ADMM' ) )
set(gca, 'FontSize', 12)
set(gca,'linewidth',2);


figure; 
set(gca,'fontsize',15)
p= plot(output_p1.transmission(:),output_p1.MSE_error(:),strcat('b-'),'Marker',Markers{3},'MarkerIndices',1:intVal2:length(output_p1.MSE_error),'MarkerSize',10,'linewidth',1.1);
hold on; plot(output_p1.transmission(:),output_p1.Z_Con_error(:),strcat('b-'),'Marker',Markers{6},'MarkerIndices',1:intVal2:length(output_p1.Z_Con_error(:)),'MarkerSize',10,'linewidth',1.1);
hold on; plot(output_p1.transmission(:),output_p1.Z_nCon_error(:),strcat('b-'),'Marker',Markers{9},'MarkerIndices',1:intVal2:length(output_p1.Z_nCon_error(:)),'MarkerSize',10,'linewidth',1.1);
hold on; plot(output_p2.transmission(:),output_p2.MSE_error(:),strcat('r-'),'Marker',Markers{3},'MarkerIndices',1:intVal2:length(output_p2.MSE_error(:)),'MarkerSize',10,'linewidth',1.1); 
hold on; plot(output_p2.transmission(:),output_p2.Z_Con_error(:),strcat('r-'),'Marker',Markers{6},'MarkerIndices',1:intVal2:length(output_p2.Z_Con_error(:)),'MarkerSize',10,'linewidth',1.1);
hold on; plot(output_p2.transmission(:),output_p2.Z_nCon_error(:),strcat('r-'),'Marker',Markers{9},'MarkerIndices',1:intVal2:length(output_p2.Z_nCon_error(:)),'MarkerSize',10,'linewidth',1.1);
ylim([error_th,Pmax^2*10])
xlim([sum(sum(aveRes.Geograph.weight)) inf])
grid on 
legend({'Proposed: x^{(k)}-x*','Proposed: (\Pi_{H})\lambda^{(k)}-\lambda*','Proposed: (I-\Pi_{H})\lambda^{(k)}','Non-private: x^{(k)}-x*','Non-private: (\Pi_{H})\lambda^{(k)}-\lambda*','Non-private: (I-\Pi_{H})\lambda^{(k)}'},'location','northeast','FontSize',10)
set(gca, 'YScale', 'log')
xlabel ('Transmissions'); ylabel ('MSE')
title(strcat('Distributed ', name,' of PDMM' ) )
set(gca, 'FontSize', 12)
set(gca,'linewidth',2);

figure; 
set(gca,'fontsize',15)
p= plot(output_d1.transmission(:),output_d1.MSE_error(:),strcat('b-'),'Marker',Markers{3},'MarkerIndices',1:intVal3:length(output_d1.MSE_error),'MarkerSize',10,'linewidth',1.1);
hold on; plot(output_d1.transmission(:),output_d1.Z_Con_error(:),strcat('b-'),'Marker',Markers{6},'MarkerIndices',1:intVal3:length(output_d1.Z_Con_error(:)),'MarkerSize',10,'linewidth',1.1);
hold on; plot(output_d1.transmission(:),output_d1.Z_nCon_error(:),strcat('b-'),'Marker',Markers{9},'MarkerIndices',1:intVal3:length(output_d1.Z_nCon_error(:)),'MarkerSize',10,'linewidth',1.1);
hold on; plot(output_d2.transmission(:),output_d2.MSE_error(:),strcat('r-'),'Marker',Markers{3},'MarkerIndices',1:intVal3:length(output_d2.MSE_error(:)),'MarkerSize',10,'linewidth',1.1); 
hold on; plot(output_d2.transmission(:),output_d2.Z_Con_error(:),strcat('r-'),'Marker',Markers{6},'MarkerIndices',1:intVal3:length(output_d2.Z_Con_error(:)),'MarkerSize',10,'linewidth',1.1);
hold on; plot(output_d2.transmission(:),output_d2.Z_nCon_error(:),strcat('r-'),'Marker',Markers{9},'MarkerIndices',1:intVal3:length(output_d2.Z_nCon_error(:)),'MarkerSize',10,'linewidth',1.1);
ylim([error_th, inf])
grid on 
legend({'Proposed: x^{(k)}-x*','Proposed: (\Pi_{H})u^{(k)}-u*','Proposed: (I-\Pi_{H})u^{(k)}','Non-private: x^{(k)}-x*','Non-private: (\Pi_{H})u^{(k)}-u*','Non-private: (I-\Pi_{H})u^{(k)}'},'location','northeast','FontSize',10)
set(gca, 'YScale', 'log')
xlabel ('Transmissions'); ylabel ('MSE')
title(strcat('Distributed ', name,' of Dual Ascent' ) )
set(gca, 'FontSize', 12)
set(gca,'linewidth',2);
