# import modules
 
import pint as pn
from pint import vmath
 
# pint is an interval type Library for Python 3.x
# https://github.com/o108minmin/pint
# distributed under MIT license

#  // verification program for the q-Ramanujan function
# // reference
# // ENCYCLOPEDIA OF MATHEMATICS AND ITS APPLICATOINS 71
# // SPECIAL FUNCTIONS
# // author: G.E.ANDREWS, RICHARD ASKEY, RANJAN ROY
# // CAMBRIDGE, 1999
# // Page 551, Exercise 39
K=500;
x=pn.interval(0.1); 
q=pn.interval(0.1);
j=1;
# // calculate Sum[q^(n^2)*(-x)^n/(q^2;q^2)_n*(xq^2;q^2)_n,{n,0,Infinity}]
a=pn.interval(1);
b=pn.interval(1);
c=pn.interval(1);
for k in range (1,K+1):
    j*=-1
    b*=1-pow(q,2*k)
    c*=1-x*pow(q,2*k)
    a+=j*pow(q,k*k)*pow(x,k)/(b*c)

b*=1-pow(q,2*K+2);
c*=1-x*pow(q,2*K+2);
rad=abs(pow(q,(K+1)*(K+1))*pow(x,K+1)/(b*c)).sup
series=a+rad;
# // calculate (xq^2;q^2)_Infinity
q1=pn.interval(q*q);
z=pn.interval(x*q1);
f=pn.interval(1/(1-z));
g=pn.interval(1);
m=1;
for l in range (1,K+1):
    m*=-1
    g*=1-pow(q1,l)
    f+=m*pow(q1,l*(l+1)/2)/(g*(1-z*pow(q1,l)))

g*=1-pow(q1,K+1)
rad2=abs(pow(q1,(K+1)*(K+2)/2)/(g*(1-z*pow(q1,K+1)))).sup
series2=f+rad2;
a1=pn.interval(1-q1);
j1=1;
for k1 in range (1,K+1):
    j1*=-1
    a1+=j1*(1-pow(q1,2*k1+1))*pow(q1,k1*(3*k1+1)/2)

rad3=abs((1-pow(q1,2*K+3))*pow(q1,(K+1)*(3*K+4)/2)).sup
euler=a1+rad3;
qexp=euler/series2;
res=qexp*series;
print(res)
