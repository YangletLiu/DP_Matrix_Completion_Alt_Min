clc;
clear all;

data = 'movielens'; %数据文件名
path = strcat('.\data\',data,'.mat');
load(path);  
D = input;   %读取数据文件
[m,n] = size(D); %返回data数据文件里的矩阵大小

RR=30;
MM=200;NN=200;
AA=randn(MM,RR);
BB=randn(RR,NN);
f01=AA*BB;
sampleee=rand(MM,NN);

sample111=sampleee';
A00=rand(MM,RR);
tic
for h=1:10
    sss=sampleee> (h*0.1);
    f02=f01.*sss;
    f021=f02';
    for t=1:10  %循环次数
        if t==1
            AA=A00; 
        end
        for j=1:NN %分别算B的每一列
            y=f02(:,j);        
            for i=1:RR
                x1(:,i)=sampleee(:,j).*AA(:,i);%将采样矩阵K的第j列跟A0对应行的每个元素相乘
            end
            b(:,j)=x1\y; %算B的第j列
        end
        BBB=b';
        for j=1:MM %分别算A'的每一列
            y=f021(:,j);
            for i=1:RR
                x2(:,i)=sample111(:,j).*BBB(:,i);
            end
            a(:,j)=x2\y; %算A'的第j列
        end
        AA=a';BBB=BBB';
        M0=AA*BBB;          
    end
    p2(1,h)=rmse(f01,M0);
end

e=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
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
toc