%%%参考文章DP_MC_Revisited的
%%%Algorithm 1 Private Frank-Wolfe algorithm
%%%分为Global Component和Local Update
addpath('./FW_T/func');
addpath('./FW_T/PROPACK');
warning off;

data = 'movielens'; %数据文件名

rho = .75;  %采样率

path = strcat('.\data\',data,'.mat');
load(path);  %读取数据文件

D = input;

[nu ni] = size(D); %返回data数据文件里的矩阵大小
fprintf('data has been loaded: m = %d, n = %d; \n', m,n);

%initialization

delta = 10^(-100);
epsilon = 2*log(1/delta);
T = 10000;
L = maxl2norm(D,rho,nu);
beta = 10^(-2);
k = 2*rank(D);

p=zeros(1,T)

sigma = (L^2*sqrt(64*T*log(1/delta)))/epsilon;
v = zeros(1*n);
lamda = 0;


for t=1:T
    (v,lamda)=global(D,rho,L,T,delta,epsilon,beta,sigma,v,lamda,k)
end
