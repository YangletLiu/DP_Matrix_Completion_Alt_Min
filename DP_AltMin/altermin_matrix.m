data = 'movielens'; %数据文件名
path = strcat('.\data\',data,'.mat');
load(path);  
D = input;   %读取数据文件
[m,n] = size(D); %返回data数据文件里的矩阵大小

rho=.75

Omega = rand(m,n)<=rho; % support of observation
D_omega = Omega.*D; % measurements are made
PP=rand(1);
P=PP(1,1);
r=rank(D);

for t=1:1000  %循环次数
    if t==1
        temp=D_omega.*P;
        [U,S1,V1]=svds(temp,R);
    end
    V=
end
p2(1,h)=rmse(f01,M0);


% subplot(1,2,1)
% imshow(uint8(f01));title('original')
% subplot(1,2,2)
%imshow(uint8(M0))
% 20*log10(p(10))
t=[1:10];
% semilogy(t,p);
semilogy(e,p2);
title('Matrix:200*200,rank:30');
xlabel('Iterations');
ylabel('RSE in log-scale');