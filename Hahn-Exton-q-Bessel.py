# import modules

import pint as pn
from pint import vmath

# pint is an interval type Library for Python 3.x
# https://github.com/o108minmin/pint
# distributed under MIT license

K=500;
# // verfication program for the Hahn-Exton q-Bessel function
j=1;
q=pn.interval(0.1);
x=pn.interval(1.4);
nu=pn.interval(3.5); # nu is a real number

a=pn.interval(1);
b=pn.interval(1);
c=pn.interval(1);

for k in range(1,K+1):
    j*=-1
    b*=1-pow(q,k)
    c*=1-pow(q,k+nu)
    a+=j*pow(q,k*(k+1)/2)*pow(x,2*k)/(b*c)

b*=1-pow(q,K+1)
c*=1-pow(q,K+1+nu)
rad=abs(pow(-1,K+1)*pow(q,(K+1)*(K+2)/2)*pow(x,2*K+2)/(b*c)).sup;
series=a+rad;
z=pn.interval(pow(q,nu+1))
f=pn.interval(1/(1-z))
g=pn.interval(1);
m=1;
for l in range(1,K+1):
    m*=-1
    g*=1-pow(q,l)
    f+=m*pow(q,l*(l+1)/2)/(g*(1-z*pow(q,l)))

g*=1-pow(q,K+1)
rad2=abs(pow(-1,K+1)*pow(q,(K+1)*(K+2)/2)/(g*(1-z*pow(q,K+1)))).sup
denom=f+rad2;
res=pow(x,nu)*series/denom;
print(res)
