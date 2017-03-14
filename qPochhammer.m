function res=qPochhammer(z,q)
%   // verification program for the q-Pochhammer symbol
%   // theorem by Fridrikh Israilevich Karpelevich was used
%   // reference: The Modified q-Bessel Functions and the q-Bessel-Macdonald Functions
%   // Olshanetsky, Rogov 1995
if q<=0
error('qPochhammer symbol only for 0<q<1')
end
if q>=1
error('qPochhammer symbol only for 0<q<1')
end

   K=500;
   KK=500;	
   n=100;
   j=1;
q=intval(q); %0<q<1

 if abs(z)<1
  a=intval(1/(1-z));
  b=intval(1);

  for k=1:K
    j = -1*j;
    b=b*(1-q^intval(k));
    a=a+j*q^(intval(k)*(intval(k)+1)/2)/(b*(1-intval(z)*q^intval(k)));
  end
  K=intval(K);
  b=b*(1-q^(K+1));
   intrad=abs(q^((K+1)*(K+2)/2)/(b*(1-intval(z)*q^(K+1))));
   format long
  series=a+midrad(0,sup(intrad)); 
% Leibniz criterion
e=eulerfunc(q);
res=intval(e/series);

else
%// Reference: Plancherel-Rotach asymptotics for certain basic hypergeometric series
%// Zhang, 2008
 while abs(intval(z))*(q^n)/(1-q)>=0.5
  n=n+500;
 end
rada=2*abs(intval(z))*(q^n)/(1-q);
format long
res=qPochhammerFinite(z,q,n)*midrad(0,sup(rada));

 if imag(z)~=0
  if (abs(z)<1) && (real(z)>0)
  x=real(z);
  y=imag(z);
  x=intval(x);
  y=intval(y);
  j1=1;
  j2=1;
  a1=(1-x)/((1-x)^2+y^2);	
  b1=1;
   for k=1:K
   j1=-1*j1;
   b1=b1*(1-q^k);
   a1=a1+j1*q^(k*(k+1)/2)*(1-x*q^k)/(b1*(1-x*q^k)^2+y*y*q^(2*k));
   end
 %Leibniz criterion
  K=intval(K);
  rad1=abs(q^((K+2)*(K+1)/2)*(1-x*q^(K+1))/(b1*(1-x*q^(K+1))^2+y*y*q^(2*K+2)));
  realz=a1+midrad(0,sup(rad1));
  a2=y/((1-x)^2+y^2);
  b2=1;
   for kk=1:KK
   j2=-1*j2;
   b2=b2*(1-q^kk);
   a2=a2+j2*q^(kk*(kk+1)/2)*(y*q^kk)/(b2*(1-x*q^kk)^2+y*y*q^(2*kk));
   end
  KK=intval(KK);
  rad2=abs(q^(KK*(KK+1)/2)*(y*q^KK)/(b2*(1-x*q^KK)^2+y*y*q^(2*KK)));
  imagz=a2+midrad(0,sup(rad2));
  format long
  res=realz+imagz*1i;
  else
%// Reference: Plancherel-Rotach asymptotics for certain basic hypergeometric series
%// Zhang, 2008
  z=intval(z);
   while abs(z)*(q^n)/(1-q)>=0.5
   n=n+500;
   end
  radb=2*abs(z)*(q^n)/(1-q);
  format long
  res=qPochhammerFinite(z,q,n)*midrad(0,sup(radb));
  end
 end
end
end

function res=qPochhammerFinite(z,q,n)
qp=intval(1);
 if n==0
  res=intval(1);
 else
  for k=0:n-1
  qp=qp*(1-z*q^k);
  end
 format long
 res=qp;
 end
end
