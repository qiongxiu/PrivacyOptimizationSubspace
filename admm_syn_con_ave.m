
function output= admm_syn_con_ave(Geograph,error_th,iteration_max,Pmax,rho,flag)
n = Geograph.n;        
N_degree = sum(Geograph.weight).';    % degree matrix (D)
numEdge=Geograph.m;
% Initialization
x_est = Geograph.node_val;
x_AverPrivate=x_est;
AverReal = mean(Geograph.node_val);
AverReal_rep = repmat(AverReal,n,1);
u_matrix=zeros(numEdge/2,1);
if flag==1|| flag==3
   v_matrix=Pmax*randn(numEdge,1);
else 
    v_matrix=zeros(numEdge,1);
end 
A_admm=Geograph.A_admm;
B_admm=Geograph.B_admm;
S_matrix=Geograph.S_admm;
S=S_matrix*pinv(S_matrix.'*S_matrix)*S_matrix.';
if flag==3
    v_matrix=S*v_matrix;
end 
Vr=S*v_matrix;
Vl=(eye(numEdge)-S)*v_matrix;
MSE_error = 1/n*(norm(x_est-AverReal*ones(n,1)))^2;
iteration = 1;
transmission=0;
tranCnt=0;%calculate the transimission times in each iteration
while iteration <= iteration_max     
     x_est=(x_AverPrivate-(A_admm.')*v_matrix-rho*A_admm.'*B_admm*u_matrix)./(1+rho*N_degree);
     u_matrix=diag(inv(B_admm.'*B_admm)).*(-1/rho*(B_admm.'*v_matrix)-B_admm.'*A_admm*x_est);
     v_matrix=v_matrix+rho*(A_admm*x_est+B_admm*u_matrix);
     MSE_error(iteration) =1/n*(norm(x_est-AverReal*ones(n,1)))^2;
     Vr1=S*v_matrix;
     Vl1=(eye(numEdge)-S)*v_matrix;
     V_nCon_error(iteration)=1/(n)*(norm(Vl1))^2;
     V_Con_error(iteration)=1/(n)*(norm(A_admm.'*Vr1-(-AverReal_rep(:)+Geograph.node_val)))^2;
     tranCnt=tranCnt+sum(N_degree);
     transmission(iteration)=tranCnt;
      if MSE_error(iteration)< error_th && V_Con_error(iteration)< error_th
          MSE_error(iteration)=error_th;
          V_Con_error(iteration)=error_th;
          V_struct.v_matrix=v_matrix;
          iteration=iteration_max+1;
      else
         iteration = iteration+1;
      end 
end
output.transmission=transmission;
output.x_est=x_est;
output.MSE_error=MSE_error;
output.Z_nCon_error=V_nCon_error;
output.Z_Con_error=V_Con_error;  
end



