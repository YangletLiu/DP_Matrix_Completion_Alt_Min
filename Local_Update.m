%%%Local Update 
%%%我们首先要将待补全的矩阵彻底的分一下，好实现分布式的搞
%%%我看文章中说的是yita是step size，也就是说可以任意的吗？？？？？
i = m;
v = ;%top right singular vector
lamda1 = ;%top singular value
T;%total number of iterations
t;%current iteration
L;%bound on l2-norm of P_\omega(Y_i^*)

Y= zeros(1*n);
A=omega(Y-Y1,rho);
u=(A.*v)/lamda1;
Y=pai((1-1/T)Y-K/T.*u.*v')
A=Omega(Y-Y1)
if T=t:
    r=Y;
else:
    r=A.*A';

