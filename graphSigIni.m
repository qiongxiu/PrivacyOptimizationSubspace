function Geograph=graphSigIni(Geograph,sigma,u,ni)
Geograph.u=u; %feature dimension
Geograph.ni=ni; %observation dimension of each node
n=Geograph.n;
if u>ni %%underdetermined system
   H_matrix=sigma*randn(n*ni,u);
   tem=randperm(u,0.1*u);
   x_true=zeros(u,1);
   x_true(tem)=randn(0.1*u,1);
   Geograph.x_true=x_true;
   Geograph.H_matrix=H_matrix; %regression matrix in the network
   Geograph.y_vec=H_matrix*Geograph.x_true; %observation vector in the network y=Hx+v, v is noise
else %%overdetermined system 
    H_matrix=sigma*randn(n*ni,u);
    Geograph.H_matrix=H_matrix; %regression matrix in the network
    Geograph.y_vec=H_matrix*randn(u,1); %observation vector in the network 
end 