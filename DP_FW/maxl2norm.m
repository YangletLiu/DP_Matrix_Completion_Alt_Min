function L=maxl2norm(D,rho,nu)

DD=omega(D,rho)
max = 0;
for i=1:nu
    buf=norm(DD(i,:));
    if max<buf
        max = buf;
    end
end

L=max;