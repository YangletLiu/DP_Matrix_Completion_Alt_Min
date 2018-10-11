function [Yi,AN]=Local_Update(i,v,lamda1,T,t,L,D,k)

Y1=D(i,:);
[~,n]=size(D);

if t==1
    Y= zeros(1*n);
end
A=omega(Y-Y1,rho);
u=(A.*v)/lamda1;
a = (1-1/T).*Y-k/T.*u.*v';
Y=projection(a,L);
A=Omega(Y-Y1);
AN=A.*A';
Yi=Y;

