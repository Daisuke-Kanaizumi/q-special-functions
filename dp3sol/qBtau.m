function res=qBtau(n,N,nu,q)
A=zeros(N);

    for k=1:N
    for l=1:N
    A(k,l)=mid(qBessel(q,nu,q^(n+k+2*l-3))); 
    % mid is used because INTLAB version 6 is not supporting LU
    end
    end
res=fastvdet(A);
% reference for determinant verification
% Ogita,Oishi,Ozaki,2007
% Oishi, Rump, 2002
% Fast verification of solutions of matrix equations
end
% alternative options for the interval of determinant
% Vdet, robustvdet