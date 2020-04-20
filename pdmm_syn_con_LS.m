
function output=pdmm_syn_con_LS(Geograph,error_th,iteration_max,Pmax,c,flag)
[Nn,u] = size(Geograph.H_matrix); 
ni=Geograph.ni;
n=Nn/ni;
N_degree = sum(Geograph.weight).';
numEdge=Geograph.m;
if flag==1||flag==3
   z_matrix=Pmax*randn(numEdge,u);
else 
    z_matrix=zeros(numEdge,u);
end 

H_matrix=Geograph.H_matrix;
y_vec=Geograph.y_vec;
C_matrix=Geograph.C_matrix;
x_est=zeros(n,u);
x_true = mldivide(H_matrix,y_vec);

Z_struct.z_ji_matrixIni=z_matrix;
S_matrix=Geograph.S_matrix;
C_matrix=Geograph.C_matrix;
PC_matrix=Geograph.PC_matrix;
S=S_matrix*pinv(S_matrix.'*S_matrix)*S_matrix.';
if flag==3
    z_matrix=S*z_matrix;
end 
Zr=S*z_matrix;
Zl=(eye(numEdge)-S)*z_matrix;
f_deri=zeros(n,u); %initialize the derivative of f(x)
for i=1:n
    ind_i=(i-1)*u+1:i*u;
    for j=1:u
        f_deri(i,j)=(H_matrix(ind_i,:)*x_true-y_vec(ind_i)).'*H_matrix(ind_i,j);
    end 
end 




Zfixed=-pinv([C_matrix.';PC_matrix.'])*[f_deri+c*(C_matrix.')*C_matrix*repmat(x_true,1,n).';f_deri+c*(C_matrix.')*PC_matrix*repmat(x_true,1,n).'];

MSE_error = 1/(n*u)*(norm(x_est-repmat(x_true,1,n).'))^2;
iteration = 1;
transmission=0;
tranCnt=0;
while iteration <= iteration_max 
    for i=1:n
         ind_i=(i-1)*ni+1:i*ni;
        x_est(i,:)=(pinv(H_matrix(ind_i,:).'*H_matrix(ind_i,:)+c*N_degree(i)*eye(u,u))*(H_matrix(ind_i,:).'*y_vec(ind_i)-(C_matrix(:,i).'*z_matrix).')).';
    end
       z_matrix=circshift(z_matrix,numEdge/2)+2*circshift(c*C_matrix*x_est,numEdge/2);

      MSE_error(iteration) = 1/(n*u)*(norm(x_est-repmat(x_true,1,n).'))^2;
      Zr1=S*z_matrix;
      Zl1=(eye(numEdge)-S)*z_matrix;
      Z_nCon_error(iteration)=1/(numEdge*u)*(norm(Zl1))^2;
      Z_Con_error(iteration)=1/(numEdge*u)*(norm(Zr1-Zfixed))^2;
      tranCnt=tranCnt+n;%each node only broadcast once a iteration 
      transmission(iteration)=tranCnt;
      if MSE_error(iteration)<error_th&& Z_Con_error(iteration)< error_th
          MSE_error(iteration)=error_th;
          Z_Con_error(iteration)=error_th;
          iteration=iteration_max+1;
          Z_struct.z_matrix=z_matrix;
      else
         iteration = iteration+1;
      end   
end
output.transmission=transmission;
output.x_est=x_est;
output.MSE_error=MSE_error;
output.Z_nCon_error=Z_nCon_error;
output.Z_Con_error=Z_Con_error;