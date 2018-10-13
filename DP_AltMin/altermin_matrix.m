data = 'movielens'; %数据文件名
path = strcat('.\data\',data,'.mat');
load(path);  
D = input;   %读取数据文件
[MM,NN] = size(D); %返回data数据文件里的矩阵大小

rho=.75

Omega = rand(MM,NN)<=rho; % support of observation
Omega1 = Omega';
D_omega = Omega.*D; % measurements are made
D_omega1 = D_omega';
PP=rand(1);
P=PP(1,1);
RR=rank(D);

for t=1:400  %循环次数
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
    VV=V';
    for j=1:MM %分别算A'的每一列
        y=D_omega1(:,j);
        for i=1:RR
            x2(:,i)=Omega1(:,j).*VV(:,i);
        end
        U(:,j)=x2\y; %算A'的第j列
    end
    UU=U';VV=VV';
    M0=UU*VV;  
    p2(1,t)=rmse(D,M0);        
end


% subplot(1,2,1)
% imshow(uint8(f01));title('original')
% subplot(1,2,2)
%imshow(uint8(M0))
% 20*log10(p(10))
%t=[1:10];
% semilogy(t,p);
%semilogy(e,p2);
%title('Matrix:200*200,rank:30');
%xlabel('Iterations');
%ylabel('RSE in log-scale');