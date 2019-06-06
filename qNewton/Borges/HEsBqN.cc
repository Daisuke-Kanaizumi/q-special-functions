#include <kv/Heine.hpp>

#include <kv/qBessel.hpp>
#include <cmath>
typedef kv::interval<double> itv;
typedef kv::complex< kv::interval<double> > cp;
using namespace std;
int main()
{
  cout.precision(17);
  int n=20;
  itv x,nu,q,x0,y;
  q="0.5";
  nu=1.5;
  x=3.3;
  x0=x;
  y=0;
  for(int i=1;i<=n;i++){
    x=x-(kv::Hahn_Exton(itv(x),itv(nu),itv(q))-y)*(1-q)*x0*(1+(1-q)*kv::Hahn_Exton(itv(q*x0),itv(nu),itv(q)))
      /(kv::Hahn_Exton(itv(x0),itv(nu),itv(q))-kv::Hahn_Exton(itv(q*x0),itv(nu),itv(q)));
  cout<<x<<endl;
  cout<<"value of HE inf"<<kv::Hahn_Exton(itv(x.lower()),itv(nu),itv(q))-y<<endl;
  cout<<"value of HE sup"<<kv::Hahn_Exton(itv(x.upper()),itv(nu),itv(q))-y<<endl;
  cout<<"value of HE mid"<<kv::Hahn_Exton(itv(mid(x)),itv(nu),itv(q))-y<<endl;
  }
  // A. L. Garcia, Numerical Methods for Physics, 2nd Edition, Pearson, 2015.
}

