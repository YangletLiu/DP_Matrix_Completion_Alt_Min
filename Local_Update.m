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
