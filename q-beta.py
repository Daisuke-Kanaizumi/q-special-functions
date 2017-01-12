# import modules

import pint as pn
from pint import vmath

# pint is an interval type Library for Python 3.x
# https://github.com/o108minmin/pint
# distributed under MIT license

# calculate the value of q-beta function

K=500;
j=1;
jj=1;
jjj=1;
print("input value of q")
# q is any real number
qq=input()
q=pn.interval(qq);
print("input value of z1")
zz1=input()
z1=pn.interval(zz1);
print("input value of z2")
zz2=input()
z2=pn.interval(zz2);
z3=pn.interval(z1+z2);
if abs(q)<1: # calculating q-beta function by using Jackson`s q-gamma function
    # calculate gamma_q(z1)
    a=pn.interval(1/(1-pow(q,z1)))
    b=pn.interval(1)
    for k in range (1,K+1):
        j*=-1
        b*=1-pow(q,k)
        a+=j*pow(q,k*(k+1)/2)/(b*(1-pow(q,z1+k)))
    
    b*=1-pow(q,K+1)
    rad=abs(pow(q,(K+1)*(K+2)/2)/(b*(1-pow(q,K+1+z1)))).sup
    series1=a+rad;
    c=pn.interval(pow(1-q,1-z1));
    res1=c*series1;
    # calculate gamma_q(z2)
    a1=pn.interval(1/(1-pow(q,z2)))
    b1=pn.interval(1)
    for kk in range (1,K+1):
        jj*=-1
        b1*=1-pow(q,kk)
        a1+=jj*pow(q,kk*(kk+1)/2)/(b1*(1-pow(q,z2+kk)))
    
    b1*=1-pow(q,K+1)
    radrad=abs(pow(q,(K+1)*(K+2)/2)/(b1*(1-pow(q,K+1+z2)))).sup
    series2=a1+radrad;
    c1=pn.interval(pow(1-q,1-z2));
    res2=c1*series2;
    # calculate gamma_q(z1+z2)
    a2=pn.interval(1/(1-pow(q,z3)))
    b2=pn.interval(1)
    for kkk in range (1,K+1):
        jjj*=-1
        b2*=1-pow(q,kkk)
        a2+=jjj*pow(q,kkk*(kkk+1)/2)/(b2*(1-pow(q,z3+kkk)))
    
    b2*=1-pow(q,K+1)
    radradrad=abs(pow(q,(K+1)*(K+2)/2)/(b2*(1-pow(q,K+1+z3)))).sup
    series3=a2+radradrad;
    c2=pn.interval(pow(1-q,1-z3));
    res3=c2*series3;
    res=res1*res2/res3;
    print(res)
elif abs(q)>1: # calculating q-beta function by using Moak`s q-gamma function
    # calculate gamma_q(z1)
    q1=pn.interval(1/qq);
    a3=pn.interval(1/(1-pow(q1,z1)))
    b3=pn.interval(1)
    j1=1;
    for k1 in range (1,K+1):
        j1*=-1
        b3*=1-pow(q1,k1)
        a3+=j1*pow(q1,k1*(k1+1)/2)/(b3*(1-pow(q1,z1+k1)))
    
    b3*=1-pow(q1,K+1)
    rad1=abs(pow(q1,(K+1)*(K+2)/2)/(b3*(1-pow(q1,K+1+z1))))
    series4=a3+rad1; 
    c3=pn.interval(pow(q-1,1-z1)*pow(q,z1*(z1-1)/2));
    res4=c3*series4;
    # calculate gamma_q(z2)
    a4=pn.interval(1/(1-pow(q1,z2)))
    b4=pn.interval(1)
    j2=1;
    for k2 in range (1,K+1):
        j2*=-1
        b4*=1-pow(q1,k2)
        a4+=j2*pow(q1,k2*(k2+1)/2)/(b4*(1-pow(q1,z2+k2)))
    
    b4*=1-pow(q1,K+1)
    rad2=abs(pow(q1,(K+1)*(K+2)/2)/(b4*(1-pow(q1,K+1+z2))))
    series5=a4+rad2; 
    c4=pn.interval(pow(q-1,1-z2)*pow(q,z2*(z2-1)/2));
    res5=c4*series5;
    # calculate gamma_q(z1+z2)
    a5=pn.interval(1/(1-pow(q1,z3)))
    b5=pn.interval(1)
    j3=1;
    for k3 in range (1,K+1):
        j3*=-1
        b5*=1-pow(q1,k3)
        a5+=j3*pow(q1,k3*(k3+1)/2)/(b5*(1-pow(q1,z3+k3)))
    
    b3*=1-pow(q1,K+1)
    rad3=abs(pow(q1,(K+1)*(K+2)/2)/(b5*(1-pow(q1,K+1+z3))))
    series6=a5+rad3; 
    c5=pn.interval(pow(q-1,1-z3)*pow(q,z3*(z3-1)/2));
    res6=c5*series6;
    res=res4*res5/res6;
    print(res)
else:
    print("this is beta function")
