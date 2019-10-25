#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <cmath>
typedef kv::interval<double> itv;
using namespace std;
int main()
{
  cout.precision(17);
  int K,n;
  K = 100;
  n=2;// order of Bessel function
  itv y,x,a,res;
  x=10.,
  y = 1.;
  for (int i = 2;i<=n;i++){
    itv ii;
    ii=i;
    y=y/ii;
  }
  a=pow(x*0.5,n)*y;
  int j;
  j = 1;
  for(int k=1;k<=K;k++){
    j =-j ;
    itv kk;
    kk=k;
    y = y/( n + kk )/kk ;
    a = a+(j*pow(x*0.5,n+2*k))*y ;
  }
  double rad;
  rad=(abs(pow(x*0.5,n+2*K+2)*y/(K+1)/(n+K+1))).upper();
  res = a + itv(-1.,1.)*rad ;
  cout<<res<<endl;
}
