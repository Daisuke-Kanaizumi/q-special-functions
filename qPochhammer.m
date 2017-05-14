function res=qPochhammer(z,q)
%   // verification program for the q-Pochhammer symbol (z,q)_infinity
%   // theorem by Fridrikh Israilevich Karpelevich was used
%   // reference: The Modified q-Bessel Functions and the q-Bessel-Macdonald Functions
%   // Olshanetsky, Rogov 1995, arXiv
% Abolished implementation
%    K=500;
%    j=1;
% q=intval(q); %-1<q<1
% z=intval(z); % |z|<1
% a=intval(1/(1-z)); 
% b=intval(1);
%  
%   for k=1:K
%     j = -1*j;
%     b=b*(1-q^intval(k));
%     a=a+j*q^(intval(k)*(intval(k)+1)/2)/(b*(1-z*q^intval(k)));
%   end
%   K=intval(K);
%   b=b*(1-q^(K+1));
%    intrad=abs(q^((K+1)*(K+2)/2)/(b*(1-z*q^(K+1))));
%    format long
%   series=a+midrad(0,sup(intrad)); 
% e=eulerfunc(q);
% res=e/series;

% refined implementation
%  // Reference: Plancherel-Rotach asymptotics for certain basic hypergeometric series
%       // Zhang, 2008 , Advances in Mathematics
if q<=0
error('qPochhammer symbol only for 0<q<1')
end
if q>=1
error('qPochhammer symbol only for 0<q<1')
end
n=100;
while abs(z)*q^n/(1-q)>=0.5
n=n+500;
end
radq=sup(2*abs(z)*q^n/(1-q));
radrad=midrad(1,radq);
res=qPn(z,q,n)*radrad;
end

function res=qPn(z,q,n)
% verification program for finite q Pochhammer symbol

if n==0
  res=1;
else
  qp=1;
  for k=0:n-1
    qp=qp*(1-z*q^k);
  end
  format long
  res=qp;
end

end
