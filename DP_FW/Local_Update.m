function r=local(i,v,lamda1,T,t,L,D,Y,k)

Y1=D(i,:);

if t==1
    Y= zeros(1*n);
end
A=omega(Y-Y1,rho);
u=(A.*v)/lamda1;
Y=projection((1-1/T)Y-k/T.*u.*v')
A=Omega(Y-Y1);
if t==Y:
    r=Y;
else:
    r=A.*A';
end

