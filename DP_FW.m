%%%参考文章DP_MC_Revisited的
%%%Algorithm 1 Private Frank-Wolfe algorithm
%%%分为Global Component和Local Update
addpath('./FW_T/func');
addpath('./FW_T/PROPACK');
warning off;

data = 'movielens'; %数据文件名

path = strcat('.\data\',data,'.mat');
load(path);  %读取数据文件

D = input;

[m n] = size(D); %返回data数据文件里的矩阵大小
fprintf('data has been loaded: m = %d, n = %d; \n', m,n);

%%%Global Component
delta = ;
epsilon = 2*log(1/delta);
iter = ;
L = ;
beta = ;
nu = m;
ni = n;