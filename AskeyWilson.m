function res=AskeyWilson(t1,t2,t3,t4,q)
% verification program for the Askey Wilson  integral
if abs(t1)>=1
    error('this program is unavailable')
elseif abs(t2)>=1
    error('this program is unavailable')
else
    res=2*qPochhammer(t1*t2*t3*t4,q)/(eulerfunc(q)*qPochhammer(t1*t2,q)*qPochhammer(t1*t3,q)*qPochhammer(t1*t4,q)*qPochhammer(t2*t3,q)*qPochhammer(t2*t4,q)*qPochhammer(t3*t4,q));
end
end
