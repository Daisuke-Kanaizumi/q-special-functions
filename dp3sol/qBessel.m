function res=qBessel(q,nu,x)
% // verfication program for the q-Bessel function
%   // reference
%   // Ultradiscrete limit of Bessel function type solutions of the the Painleve III equation
%   // Isojima, 2014, arXiv
%   // final goal
%   // verify the determinantial structure solution of discrete Painleve III equation
format long
K=2000;
j=1;
% q must be positive, |q|<1
% nu is a positive integer
q=intval(q);
nu=intval(nu);
x=intval(x);
b=intval(1);
c=intval(1);
for l=1:mid(nu)
    c=c*(1-q^(2*l));
end
a=intval(1/c);
for k=1:K
    j=-1*j;
    b=b*(1-q^(2*k));
    c=c*(1-q^(2*(k+nu)));
    a=a+j*(x^(2*k))*((1-q)^(2*k))/(b*c);
end
K=intval(K);
b=b*(1-q^(2*(K+1)));
c=c*(1-q^(2*(K+nu+1)));
rada=abs(((x/2)^(2*K+2))*((1-q)^(2*K+2))/(b*c));
series=a+rada;
res=((1-q)^(nu))*(x^(nu))*series;
end