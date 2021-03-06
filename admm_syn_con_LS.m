function output=admm_syn_con_LS(Geograph,error_th,iteration_max,Pmax,c,flag)
[Nn,u] = size(Geograph.H_matrix); 
ni=Geograph.ni;
n=Nn/ni;
N_degree = sum(Geograph.weight).';
numEdge=Geograph.m;
u_matrix=zeros(numEdge/2,u);
if flag==1||flag==3
   v_matrix=Pmax*randn(numEdge,u);
else 
    v_matrix=zeros(numEdge,u);
end 
V_struct.v_ji_matrixIni=v_matrix; 
A_admm=Geograph.A_admm;
B_admm=Geograph.B_admm;
S_matrix=Geograph.S_admm;
S=S_matrix*pinv(S_matrix.'*S_matrix)*S_matrix.';
if flag==3
    v_matrix=S*v_matrix;
end 
Vr=S*v_matrix;
Vl=(eye(numEdge)-S)*v_matrix;
H_matrix=Geograph.H_matrix;
y_vec=Geograph.y_vec;
x_true = mldivide(H_matrix,y_vec);
x_est=zeros(n,u);
for i=1:n
    ind_i=(i-1)*ni+1:i*ni;
    for j=1:u
        f_deri(i,j)=(H_matrix(ind_i,:)*x_true-y_vec(ind_i)).'*H_matrix(ind_i,j);
    end 
end 


MSE_error = 1/(n*u)*(norm(x_est-repmat(x_true,1,n).'))^2;
iteration = 1;
transmission=0;
tranCnt=0;
while iteration <= iteration_max 
    for i=1:n
      ind_i=(i-1)*ni+1:i*ni;
      x_est(i,:)=(pinv(H_matrix(ind_i,:).'*H_matrix(ind_i,:)+c*N_degree(i)*eye(u,u))*(H_matrix(ind_i,:).'*y_vec(ind_i)-(c*(A_admm(:,i).')*B_admm*u_matrix+A_admm(:,i).'*v_matrix).')).';
    end 
    u_matrix=diag(inv(B_admm.'*B_admm)).*(-1/c*(B_admm.'*v_matrix)-B_admm.'*A_admm*x_est);
    v_matrix=v_matrix+c*(A_admm*x_est+B_admm*u_matrix);
    MSE_error(iteration) = 1/(n*u)*(norm(x_est-repmat(x_true,1,n).'))^2;
    Vr1=S*v_matrix;
      Vl1=(eye(numEdge)-S)*v_matrix;
      V_nCon_error(iteration)=1/(numEdge*u)*(norm(Vl1))^2;
      V_Con_error(iteration)=1/(numEdge*u)*(norm(A_admm.'*Vr1-(f_deri)))^2;
         tranCnt=tranCnt+numEdge;
     transmission(iteration)=tranCnt;
      if MSE_error(iteration)<error_th && V_Con_error(iteration)< error_th
          MSE_error(iteration)=error_th;
          V_Con_error(iteration)=error_th;
          iteration=iteration_max+1;
          V_struct.v_matrix=v_matrix;
      else
         iteration = iteration+1;
      end   
end
output.transmission=transmission;
output.x_est=x_est;
output.MSE_error=MSE_error;
output.Z_nCon_error=V_nCon_error;
output.Z_Con_error=V_Con_error;