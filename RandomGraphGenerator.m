%%%%routing of multipoint connections%%%%
function Geograph = RandomGraphGenerator(Room_size,TotalNodeNum,sigma,sign)
radius=sqrt(2*log(TotalNodeNum)/TotalNodeNum);
node_x=Room_size(1).*abs(rand(TotalNodeNum,1));
node_y=Room_size(1).*abs(rand(TotalNodeNum,1));
A=ones(TotalNodeNum);
B=ones(TotalNodeNum);
A=A-diag(diag(B));
for i=1:TotalNodeNum
    for j=i+1:TotalNodeNum
        d=sqrt(power(node_x(i)-node_x(j),2)+power(node_y(i)-node_y(j),2)); 
        if radius<d
            A(i,j)=0;
            A(j,i)=0;
        end  
    end 
end
Geograph.n=TotalNodeNum;
Geograph.graph = graph(A~=0);
Geograph.weight=A;%adjacent matrix n*n
%initial state value of average consensus 
Geograph.node_val=randn(TotalNodeNum,1);

Inci_matrix = full(incidence(Geograph.graph ))';%incidence matrix numEdge/2 * n
Geograph.Inci_matrix=Inci_matrix;
numEdge=sum(sum(A));%number of directed edges, two times of the number of undirected edges 
Geograph.m=numEdge;

%PDMM related 
Geograph.C_matrix=[-1*(Inci_matrix<0);Inci_matrix>0]; %C matrix in PDMM
Geograph.PC_matrix = circshift(Geograph.C_matrix,numEdge/2); %PC matrix in PDMM
Geograph.S_matrix=[Geograph.C_matrix,Geograph.PC_matrix]; %incidence matrix of bipartite graph in PDMM

%ADMM related 
Geograph.A_admm=[abs(Inci_matrix<0);abs(Inci_matrix>0)];
B_matrix=[-1*eye(numEdge/2,numEdge/2);-1*eye(numEdge/2,numEdge/2)];
Geograph.B_admm=B_matrix;
Geograph.S_admm=[Geograph.A_admm,Geograph.B_admm];%incidence matrix of bipartite graph in ADMM

if rank(Inci_matrix)==TotalNodeNum-1 
      sign=1;
else 
    print('Graph is not connected, please generate again')
end 
Geograph.sign=sign;%sign=1 means the graph is connected








