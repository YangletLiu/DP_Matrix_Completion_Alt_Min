%%%Global Component
%initialization
delta = 10^(-100);
epsilon = 2*log(1/delta);
T = ;
L = maxl2norm(D,rho,nu);
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