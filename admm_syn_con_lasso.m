
function output=admm_syn_con_lasso(Geograph,error_th,iteration_max,Pmax,c,lamda,flag)
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
v_ji_matrixIni=v_matrix; 
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
x_true =Geograph.x_true;
x_est=zeros(n,u);
x_dif=x_true;
x_dif(find(x_true>0),:)=1;
x_dif(find(x_true<0),:)=-1;
for i=1:n
    ind_i=(i-1)*ni+1:i*ni;
    for j=1:u
        f_deri(i,j)=(H_matrix(ind_i,:)*x_true-y_vec(ind_i)).'*H_matrix(ind_i,j);
    end 
    f_deri(i,:)=f_deri(i,:)+x_dif';
end 
%%save the matrix inverse result
invH=zeros(n,u,u);
for i=1:n
   ind_i=(i-1)*ni+1:i*ni;
    invH(i,:,:)=pinv(H_matrix(ind_i,:).'*H_matrix(ind_i,:)+c*N_degree(i)*eye(u,u));
end 

tem=A_admm*x_est+B_admm*u_matrix;
MSE_error = (norm(repmat(y_vec,1,n)-H_matrix*x_est.')+lamda*norm(x_est,1)+sum(diag(v_matrix*(tem)'))+0.5*c*sum(sum(tem.^2,2).^(1/2)));
iteration = 1;
transmission=0;
tranCnt=0;
while iteration <= iteration_max 
    for i=1:n
      ind_i=(i-1)*ni+1:i*ni;
      x_est(i,:)=(reshape(invH(i,:,:),u,u)*wthresh(H_matrix(ind_i,:).'*y_vec(ind_i)-(c*(A_admm(:,i).')*B_admm*u_matrix+A_admm(:,i).'*v_matrix).','s',lamda)).';
    end 
    u_matrix=diag(inv(B_admm.'*B_admm)).*(-1/c*(B_admm.'*v_matrix)-B_admm.'*A_admm*x_est);
    v_matrix=v_matrix+c*(A_admm*x_est+B_admm*u_matrix);
     MSE_error(iteration) = 1/(n*u)*(norm(x_est-repmat(x_true,1,n).'))^2;
    Vr1=S*v_matrix;
      Vl1=(eye(numEdge)-S)*v_matrix;
      V_nCon_error(iteration)=1/(n*u)*(norm(Vl1))^2;
      V_Con_error(iteration)=1/(n*u)*(norm(A_admm.'*Vr1-(f_deri)))^2;
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
end 

