% This is test function for 'Privacy-preserving distributed optimization via subspace
% perturbation: a general framework'
% written by Qiongxiu Li, qili@create.aau.dk, in Aalborg, Denmark on 01.02.2020.


%%%test function for three applications: avergae, least square, lasso
close all;clear all;clc;
%%%Global parameter setting 
TotalNodeNum=20;%Network size 
c=0.1;%penalty parameter of ADMM and PDMM
Room_size = [1,1]; 
sigma=1;%variance of private data 
 
%%%Random geometric graph generaation
sign=0; 
Geograph= RandomGraphGenerator(Room_size,TotalNodeNum,sigma,sign);
Geograph_ave=Geograph;

flag=0;%0 means zero initalization, 1 means random initialization, 3 means nonzero initializatin within converging subspace 
% %%%distributed average consensus application 
fprintf('Generating results for distributed average consensus\n')
error_th = 1e-10;
iteration_max = 1e3;
Pmax=1000*sigma; %noise variance 
ave_Res.Geograph=Geograph_ave;
ave_Res.intVal=[12,7,5];
ave_Res.output_a1=admm_syn_con_ave(Geograph_ave,error_th,iteration_max,Pmax,c,1);
ave_Res.output_a2=admm_syn_con_ave(Geograph_ave,error_th,iteration_max,Pmax,c,3);
ave_Res.output_p1=pdmm_syn_con_ave(Geograph_ave,error_th,iteration_max,Pmax,c,1);
ave_Res.output_p2=pdmm_syn_con_ave(Geograph_ave,error_th,iteration_max,Pmax,c,3);
ave_Res.output_d1=dual_syn_con_ave(Geograph_ave,error_th,iteration_max,Pmax,c,1);
ave_Res.output_d2=dual_syn_con_ave(Geograph_ave,error_th,iteration_max,Pmax,c,3);
subPert_converplot(error_th,ave_Res,Pmax,' Average Consensus');

ave_Res.admm(1)=ave_Res.output_a1;
ave_Res.pdmm(1)=ave_Res.output_p1;
ave_Res.dual(1)=ave_Res.output_d1;
ave_Res.intVal=[10,5,5];
ave_Res.dp(1) = Synchronous_dp( Geograph_ave, error_th,iteration_max,0);
for i=1:2
    Pmax=Pmax/10;
    ave_Res.admm(i+1)=admm_syn_con_ave(Geograph_ave,error_th,iteration_max,Pmax,c,1);
    ave_Res.pdmm(i+1)=pdmm_syn_con_ave(Geograph_ave,error_th,iteration_max,Pmax,c,1);
    ave_Res.dual(i+1)=dual_syn_con_ave(Geograph_ave,error_th,iteration_max,Pmax,c,1);
    ave_Res.dp(i+1) = Synchronous_dp( Geograph_ave, error_th,iteration_max,Pmax);
end 
subPert_privplot(error_th,ave_Res,' Average Consensus');


% % distributed least square application without regularization 
fprintf('Generating results for distributed least squares\n')
u=3;ni=5;  %overdetermined system
Geograph_LS=graphSigIni(Geograph,sigma,u,ni);
Pmax=1000*sigma;
LS_Res.Geograph=Geograph_LS;
error_th = 1e-10;
iteration_max = 1e4;
LS_Res.intVal=[200,80,1000];
LS_Res.output_a1=admm_syn_con_LS(Geograph_LS,error_th,iteration_max,Pmax,c,1);
LS_Res.output_a2=admm_syn_con_LS(Geograph_LS,error_th,iteration_max,Pmax,c,3);
LS_Res.output_p1=pdmm_syn_con_LS(Geograph_LS,error_th,iteration_max,Pmax,c,1);
LS_Res.output_p2=pdmm_syn_con_LS(Geograph_LS,error_th,iteration_max,Pmax,c,3);
LS_Res.output_d1=dual_syn_con_LS(Geograph_LS,error_th,iteration_max,Pmax,c,1);
LS_Res.output_d2=dual_syn_con_LS(Geograph_LS,error_th,iteration_max,Pmax,c,3);
subPert_converplot(error_th,LS_Res,Pmax,' Least Squares');
% 
LS_Res.admm(1)=LS_Res.output_a1;
LS_Res.pdmm(1)=LS_Res.output_p1;
LS_Res.dual(1)=LS_Res.output_d1;
for i=1:2
    Pmax=Pmax/10;
    LS_Res.admm(i+1)=admm_syn_con_LS(Geograph_LS,error_th,iteration_max,Pmax,c,1);
    LS_Res.pdmm(i+1)=pdmm_syn_con_LS(Geograph_LS,error_th,iteration_max,Pmax,c,1);
    LS_Res.dual(i+1)=dual_syn_con_LS(Geograph_LS,error_th,iteration_max,Pmax,c,1);
    
end 
subPert_privplot(error_th,LS_Res,' Least Squares');



%%%distributed least square with L1 regularization, Lasso
fprintf('Generating results for distributed lasso\n')
u=100;ni=3; %underdetermined system
Geograph_lasso=graphSigIni(Geograph,sigma,u,ni);
Pmax=1000*sigma;
error_th = 1e-3;
c=0.2;%penalty parameter of ADMM and PDMM
iteration_max = 1e3;
lamda=0.4;%%sparsity controller
Lasso_Res.intVal=[120,180,120];
Lasso_Res.Geograph=Geograph_lasso;
Lasso_Res.output_a1=admm_syn_con_lasso(Geograph_lasso,error_th,iteration_max,Pmax,c,lamda,1);
Lasso_Res.output_a2 =admm_syn_con_lasso(Geograph_lasso,error_th,iteration_max,Pmax,c,lamda,3);
Lasso_Res.output_d1=dual_syn_con_lasso(Geograph_lasso,error_th,iteration_max,Pmax,c,lamda,1);
Lasso_Res.output_d2=dual_syn_con_lasso(Geograph_lasso,error_th,iteration_max,Pmax,c,lamda,3);
iteration_max = 2e3;
Lasso_Res.output_p1=pdmm_syn_con_lasso(Geograph_lasso,error_th,iteration_max,Pmax,c,lamda,1);
Lasso_Res.output_p2=pdmm_syn_con_lasso(Geograph_lasso,error_th,iteration_max,Pmax,c,lamda,3);
subPert_converplot(error_th,Lasso_Res,Pmax,' LASSO');
% Lasso_Res.admm(1)=Lasso_Res.output_a1;
% Lasso_Res.pdmm(1)=Lasso_Res.output_p1;
% Lasso_Res.dual(1)=Lasso_Res.output_d1;
% for i=1:2
%     Pmax=Pmax/10;
%     Lasso_Res.admm(i+1)=admm_syn_con_lasso(Geograph_lasso,error_th,iteration_max,Pmax,c,lamda,1);
%     Lasso_Res.pdmm(i+1)=pdmm_syn_con_lasso(Geograph_lasso,error_th,iteration_max,Pmax,c,lamda,1);
%     Lasso_Res.dual(i+1)=dual_syn_con_lasso(Geograph_lasso,error_th,iteration_max,Pmax,c,lamda,1);
% end 
% subPert_privplot(error_th,Lasso_Res,' LASSO');


% %%%comparison with differential privacy for distributed average consensus 
fprintf('Generating results for comarison with differential privacy in distributed average consensus\n')
error_th = 1e-5;
iteration_max = 2e2;
Pmax=sigma;
ave_Res.Geograph=Geograph_ave;
ave_Res.intVal=[12,7,5];
ave_Res.pdmm(1)=pdmm_syn_con_ave(Geograph_ave,error_th,iteration_max,0,c,1);
ave_Res.admm(1)=admm_syn_con_ave(Geograph_ave,error_th,iteration_max,0,c,1);
ave_Res.dp(1) = Synchronous_dp( Geograph_ave, error_th,iteration_max,0);
for i=1:2
    Pmax=Pmax*10;
    ave_Res.pdmm(i+1)=pdmm_syn_con_ave(Geograph_ave,error_th,iteration_max,Pmax,c,1);
    ave_Res.admm(i+1)=admm_syn_con_ave(Geograph_ave,error_th,iteration_max,Pmax,c,1);
    ave_Res.dp(i+1) = Synchronous_dp( Geograph_ave, error_th,iteration_max,Pmax);
end 
Markers = {'o','x','s','v','d','^','>','<'};
figure; 
set(gca,'fontsize',12)
for i=1:3
    plot( ave_Res.pdmm(i).transmission(:),ave_Res.pdmm(i).MSE_error(:),strcat('b-'),'Marker',Markers{i},'MarkerIndices',1:5:length(ave_Res.pdmm(i).MSE_error),'MarkerSize',10,'linewidth',1.1);
     hold on 
     plot( ave_Res.admm(i).transmission(:),ave_Res.admm(i).MSE_error(:),strcat('g-'),'Marker',Markers{i},'MarkerIndices',1:10:length(ave_Res.admm(i).MSE_error),'MarkerSize',10,'linewidth',1.1);
    hold on
     plot( ave_Res.dp(i).transmission(:),ave_Res.dp(i).MSE_error(:),strcat('r-'),'Marker',Markers{i},'MarkerIndices',1:20:length(ave_Res.dp(i).MSE_error),'MarkerSize',10,'linewidth',1.1);
    hold on 
end
set(gca,'fontsize',15)
ylim([error_th 1e5])
yticks([error_th 1e1 1e5])
xlim([0 3e4])
grid on 
legend({'p-PDMM: $10^4$','p-ADMM: $10^4$','DP: $10^4$','p-PDMM: $10^2$','p-ADMM: $10^2$','DP: $10^2$','p-PDMM: 0','p-ADMM: 0','DP: $0$'},'location','best','FontSize',15,'Interpreter','Latex')
set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
xlabel ('Transmissions'); ylabel ('||x^{(k)}-x*||^{2}')
set(gca, 'FontSize', 15)
set(gca,'linewidth',1.1);




