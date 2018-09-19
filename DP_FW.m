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

%%%Local Update 
%%%我们首先要将待补全的矩阵彻底的分一下，好实现分布式的搞
%%%我看文章中说的是yita是step size，也就是说可以任意的吗？？？？？
i = m;
v = ;
lamda1 = ;
T;
t;
L;

%根据rho来确定采样Omega，实现投影操作
if rho == 1
    
    fprintf('RPCA with full obseravation; \n');
    obs = D; Omega = ones(m,n);
    
else
    
    fprintf('RPCA with partial obseravation: ');
    Omega = rand(m,n)<=rho; % support of observation
    obs = Omega.*D; % measurements are made
    fprintf('observations are generated; \n');
    
end



%%%Global Component
%initialization
delta = ;
epsilon = 2*log(1/delta);
T = ;
L = ;
beta = ;
m = nu;
n = ni;

sigma = (L^2*sqrt(64*iter*log(1/delta)))/epsilon;
v = zeros(1*n);
lamda = 0;
for t = 1:T
    W = zeros(n*n);
    lamda1 = lamda + sqrt(sigma*log(n/delta))*n^(1/4);
    for i = 1:m
       W = W + %local(i,v,lamda,T,t,L)% 
    end
    W = W + normrnd(0,sigma^2);
    (v,lamda) = eig(W);
end