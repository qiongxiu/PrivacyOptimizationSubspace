
function output = pdmm_syn_con_ave(Geograph, error_th, iteration_max,Pmax,rho,flag)
n = Geograph.n;    
N_degree = sum(Geograph.weight).';    % degree matrix (D)
% Initialization
x_AverEst = Geograph.node_val;
x_AverPrivate=x_AverEst;
AverReal = mean(Geograph.node_val);
c= rho;
A=Geograph.weight;
numEdge=Geograph.m;
if flag==1|| flag==3
   z_ji_vec=Pmax*randn(numEdge,1);
else 
    z_ji_vec=zeros(numEdge,1);
end 
Z_struct.z_ji_matrixIni=z_ji_vec;
S_matrix=Geograph.S_matrix;
C_matrix=Geograph.C_matrix;
PC_matrix=Geograph.PC_matrix;
S=S_matrix*pinv(S_matrix.'*S_matrix)*S_matrix.';
if flag==3
    z_ji_vec=S*z_ji_vec;
end 

Zr=S*z_ji_vec;
Zl=(eye(numEdge)-S)*z_ji_vec;
AverReal_rep = repmat(AverReal,n,1);
Zfixed=-pinv([C_matrix.';PC_matrix.'])*[AverReal_rep(:)-Geograph.node_val+c*(C_matrix.')*C_matrix*AverReal_rep(:);AverReal_rep(:)-Geograph.node_val+c*(C_matrix.')*PC_matrix*AverReal_rep(:)];

MSE_error = 1/(n)*(norm(x_AverEst(:)-AverReal_rep(:)))^2;
iteration = 1;
transmission=0;
tranCnt=sum(N_degree);
% tranCnt=sum(sum(Geograph.weight));%calculate the transimission times in each iteration

while iteration <= iteration_max     
     x_est=(x_AverPrivate-(C_matrix.')*z_ji_vec)./diag(eye(n)+c*(C_matrix.')*C_matrix);
     z_ji_vec=2*circshift(c*C_matrix*x_est,numEdge/2)+circshift(z_ji_vec,numEdge/2);
      x_AverEst=x_est;
      MSE_error(iteration) = 1/(n)*(norm(x_AverEst(:)-AverReal_rep(:)))^2;
      Zr1=S*z_ji_vec;
      Zl1=(eye(numEdge)-S)*z_ji_vec;
      Z_nCon_error(iteration)=1/(numEdge)*(norm(Zl1))^2;
      Z_Con_error(iteration)=1/(numEdge)*(norm(Zr1-Zfixed))^2;
      tranCnt=tranCnt+n;%each node only broadcast once a iteration 
      transmission(iteration)=tranCnt;
      if MSE_error(iteration)<error_th && Z_Con_error(iteration)< error_th
          MSE_error(iteration)=error_th;
          Z_Con_error(iteration)=error_th;
          iteration=iteration_max+1;
          Z_struct.z_ji_vec=z_ji_vec;
      else
         iteration = iteration+1;
      end  
end
output.transmission=transmission;
output.x_est=x_est;
output.MSE_error=MSE_error;
output.Z_nCon_error=Z_nCon_error;
output.Z_Con_error=Z_Con_error;
