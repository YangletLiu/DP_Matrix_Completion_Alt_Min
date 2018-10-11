function r=omega(D,rho)

[m,n]=size(D);
if rho == 1
    
    %fprintf('RPCA with full obseravation; \n');
    Omega = ones(m,n);
    r = Omega.*D; % measurements are made
     %fprintf('observations are generated; \n');
    
else
    
    %fprintf('RPCA with partial obseravation: ');
    Omega = rand(m,n)<=rho; % support of observation
    r = Omega.*D; % measurements are made
    %fprintf('observations are generated; \n');

end