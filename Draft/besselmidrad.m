function res=besselmidrad(n,x)
x = intval(x);
n = intval(n);
%K=ceil(sup((1/2)*(sqrt(n^2+x^2)-n))); limited to double
K = 100;
y = intval(1);
for i = 2:mid(n)
    y = y./ i;
end
%factorial n
a=intval(((x./2).^n).*y);
%y = intval(0.5);
format long
j = 1;
for k=1:K
    j = -1.*j;
    y = y .* 1./(n + k)./(intval(k));
%    a=a+(((-1)^k)*(x/2)^(n+2*intval(k)))/(factorial(k)*factorial(n+k));
    a = a+(j.*(x./2).^(n+2.*intval(k))).*y;
end
%mida=a
%disp('middle')
%disp(mida)
K = intval(K);
rada=abs(intval((((-1).^(K+1)).*(x./2).^(n+2.*K+2)).*y./(K+1)./(n+K+1)));

res = a + midrad(0,sup(rada));

end
% next goal: make this program work for any x
% intervals must be controlled