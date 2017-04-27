function res=qPn(z,q,n)
% verification program for finite q Pochhammer symbol
qp=1;
if n==0
  res=1;
else
for k=0:n-1
    qp=qp*(1-z*q^k);
end
format long
res=qp;
end

