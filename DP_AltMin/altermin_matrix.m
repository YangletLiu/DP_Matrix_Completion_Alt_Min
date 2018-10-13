data = 'movielens'; %数据文件名
path = strcat('.\data\',data,'.mat');
load(path);  
D = input;   %读取数据文件
[MM,NN] = size(D); %返回data数据文件里的矩阵大小

rho=.75;

Omega = rand(MM,NN)<=rho; % support of observation
D_omega = Omega.*D; % measurements are made
PP=rand(1);
P=PP(1,1);
RR=rank(D);
T=10;

for t=1:T  %循环次数
    if t==1
        temp=D_omega.*P;
        [U,S1,V1]=svds(temp,RR);
    end
    for j=1:NN %分别算V的每一列
        y=D_omega(:,j);        
        for i=1:RR
            x1(:,i)=Omega(:,j).*U(:,i);%将采样矩阵K的第j列跟A0对应行的每个元素相乘
        end
        V(:,j)=x1\y; %算B的第j列
    end
    for j=1:MM %分别算U的每一行
        y=D_omega(j,:);
        for i=1:RR
            x2(i,:)=Omega(j,:).*V(i,:);
        end
        U(j,:)=y/x2; %算A'的第j列
    end
    M0=U*V;  
    p2(1,t)=rmse(D,M0);        
end

q=[1:T];
pic=semilogy(q,p2);
imwrite(pic,picture .fig');
I=imread(picture.fig);
imwrite(picture1.bmp);
imwrite(picture1.png);
