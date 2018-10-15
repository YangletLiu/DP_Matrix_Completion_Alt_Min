data = 'movielens'; %数据文件名
path = strcat('.\data\',data,'.mat');
load(path);  
D = input;   %读取数据文件
[MM,NN] = size(D); %返回data数据文件里的矩阵大小

rho=.75;

Omega = rand(MM,NN)<=rho; % support of observation
D_omega = Omega.*D; % measurements are made
D_omega1= D_omega';
Omega1 =Omega';
% PP=rand(1);
% P=PP(1,1);
RR=rank(D);
T=10;


I = maxl2norm(D,Omega);
delta=2.2251e-308;
temp=sqrt(2*log(2/(delta)))/(2*log(1/(delta)));
sigma = 2*I*sqrt(2*log(2/(delta)))/(2*log(1/(delta)));

for t=1:T  %循环次数
    if t==1
        U = randi([0,5],MM,RR);
%         temp=D_omega.*P;
%         [U,S1,V1]=svds(temp,RR);
    end
    for j=1:NN %分别算V的每一列
        y=D_omega(:,j);        
        for i=1:RR
            x1(:,i)=Omega(:,j).*U(:,i);%将采样矩阵K的第j列跟U0对应行的每个元素相乘
        end
        v(:,j)=pinv(x1)*y; %算v的第j列
    end
    V=v';
    for j=1:MM %分别算U的每一行
        y=D_omega1(:,j);
        for i=1:RR
            x2(:,i)=Omega1(:,j).*V(:,i);
        end
        u(:,j)=pinv(x2)*y; %算A'的第j列
    end
    U=u';V=V';
    U=U+normrnd(0,sigma^2);
    M0=U*V;  
    p(t)=norm(D-M0,'fro')/norm(D(:));
    p2(t)=rmse(D,M0);        
end


% 20*log10(p(10))
q=[1:T];
pic=semilogy(q,p2);
