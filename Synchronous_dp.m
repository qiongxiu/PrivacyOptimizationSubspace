%%%%%% paper  ''Differentially private average consensus: Obstructions, trade-offs, and
%%%%%%optimal algorithm design'' by Nozari.et.al 

function output = Synchronous_dp( Geograph,  error_th, iteration_max,sigma) 
[M,numb_aver] = size(Geograph.node_val);        
x_AverEst = Geograph.node_val;
AverReal = mean(Geograph.node_val);
AverReal_matrix = ones(M,1)*AverReal;
MSE_error = 1/(M*numb_aver)*(norm(x_AverEst(:)-AverReal_matrix(:)))^2;
iteration = 1;
tranCnt=0;%calculate the transimission times in each iteration
dmax=max(sum(Geograph.weight));
E=Geograph.weight;
a=1/dmax;
h=a/2; %h should be smaller than a 
L=diag(sum(E))-E;%graph laplacian matrix 
% s=2*rand(M,1);
s=ones(M,1);
S=diag(s);
q_vec=abs(s-1)+(1-abs(s-1)).*rand(M,1);
thetaVec=x_AverEst;
sum_vec=zeros(M,1);
numEdge=sum(sum(Geograph.weight));
while iteration <= iteration_max
    %%% privacy concerned sycronous averaing using differential privacy 
       r_vec=sigma*randl(M,1).*power(q_vec,iteration);
       x_vec=thetaVec+r_vec;%update transmitted message with laplacian nosie with varaince sigma 
%        thetaVec=thetaVec-h*L*x_vec;
       thetaVec=thetaVec-h*L*x_vec+S*r_vec;
       sum_vec=sum_vec+r_vec;
       x_AverEst=thetaVec;
       MSE_error(iteration) = 1/(M*numb_aver)*(norm(x_AverEst(:)-AverReal_matrix(:)))^2;
    tranCnt=tranCnt+numEdge;
    transmission(iteration)=tranCnt;
    if MSE_error(iteration)<error_th
          MSE_error(iteration)=error_th;
          iteration=iteration_max+1;
          iteration = iteration+1;
    else
         iteration = iteration+1;
    end 
    
end
output.transmission=transmission;
output.x_est=x_AverEst;
output.MSE_error=MSE_error;




