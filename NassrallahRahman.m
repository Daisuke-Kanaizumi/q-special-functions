function res=NassrallahRahman(t0,t1,t2,t3,t4,q)
% verification program for the Nassrallah Rahman integral
if abs(t0)>=1 || abs(t1)>=1 || abs(t2)>=1 || abs(t3)>=1 || abs(t4)>=1
    error('this program is unavailable')
else
    A=t0*t1*t2*t3*t4;
   num=2*qPochhammer(A*t0^(-1),q)*qPochhammer(A*t1^(-1),q)*qPochhammer(A*t2^(-1),q)*qPochhammer(A*t3^(-1),q)*qPochhammer(A*t4^(-1),q);
   denom=eulerfunc(q)*qPochhammer(t0*t1,q)*qPochhammer(t0*t2,q)*qPochhammer(t0*t3,q)*qPochhammer(t0*t4,q)*qPochhammer(t1*t2,q)*qPochhammer(t1*t3,q)*qPochhammer(t1*t4,q)*qPochhammer(t2*t3,q)*qPochhammer(t2*t4,q)*qPochhammer(t3*t4,q);
res=num/denom;
end
end
% reference
% Hyperbolic Beta Integrals, Jasper V. Stokman