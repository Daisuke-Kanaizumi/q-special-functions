function res=eulerfunc(q)
% verification program for Euler function
q=intval(q);
% remark : 0<|q|<1
K=500;
j=1;
a=intval(1-q);
for k=1:K
    j=-1*j;
    a=a+j.*(1-q^(2*intval(k)+1))*q^(intval(k)*(3*intval(k)+1)./2);
end
K=intval(K);
rada=abs(intval((-1).^(K+1).*(1-q^(2*K+3))*q^((K+1)*(3*K+4)/2)));
format long
res=a+midrad(0,sup(rada));
end