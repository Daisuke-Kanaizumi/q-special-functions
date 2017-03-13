function res=dp3sol(n,N,nu,q)
% Verification program for the particular solution for dPIII eq
% The solution has a determinant structure
% Components are q-Bessel functions
%    reference
%    Ultradiscrete limit of Bessel function type solutions of the the Painleve III equation
%    Isojima, 2014, arXiv
a=qBtau(n+1,N+1,nu,q);
a1=qBtau(n,N,nu+1,q);
b=qBtau(n,N+1,nu,q);
b1=qBtau(n+1,N,nu+1,q);
res=a*a1/(b*b1)-q^(nu+N);
end