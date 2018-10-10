function (v1,lamda2)=global(D,rho,L,T,delta,epsilon,beta,sigma,v,lamda)
%%%Global Component


W = zeros(n*n);
lamda1 = lamda + sqrt(sigma*log(n/beta))*n^(1/4);
for i = 1:m
    W = W + local(i,v,lamda1,T,t,L,D,Y) 
end
W = W + normrnd(0,sigma^2);
(v,lamda) = eig(W);
