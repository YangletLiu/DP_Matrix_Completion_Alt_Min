function [Yi,AN]=Local_Update(i,v,lamda1,T,t,L,D,k,Omega,Yi)

Y1=D(i,:);
Omegai=Omega(i,:);
[~,n]=size(D);

if t==1
    Y= zeros(1,n);
else
    Y=Yi;
end
A=omega(Y-Y1,Omegai);
u=(A.*v)/lamda1;
a = (1-1/T)*Y-k/T*u*v';
Y=projection(a,L,Omegai);
A=omega(Y-Y1,Omegai);
AN=A.*A';
Yi=Y;

