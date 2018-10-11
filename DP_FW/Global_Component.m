data = 'movielens'; %数据文件名
path = strcat('.\data\',data,'.mat');
load(path);  
D = input;   %读取数据文件
[m,n] = size(D); %返回data数据文件里的矩阵大小
%fprintf('data has been loaded: m = %d, n = %d; \n', m,n);

%初始化global部分所需的参数
rho = .75;  %采样率
delta = 10^(-100);
epsilon = 2*log(1/delta);
T = 10000;
L = maxl2norm(D,rho,m);
beta = 10^(-2);
k = 2*rank(D);

p=zeros(1,T);

sigma = (L^2*sqrt(64*T*log(1/delta)))/epsilon;
v = zeros(1*n);
lamda = 0;


for t=1:T
    W = zeros(n*n);
    lamda1 = lamda + sqrt(sigma*log(n/beta))*n^(1/4);
    for i = 1:m
        W = W + local(i,v,lamda1,T,t,L,D,Y) ;
    end
    W = W + normrnd(0,sigma^2);
    [V, D] = eig(W);
    lambda = wrev(diag(D));
    V = fliplr(V);
end
