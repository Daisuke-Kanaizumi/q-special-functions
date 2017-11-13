#include <kv/Heine.hpp>
#include <kv/qAiry.hpp>
#include <kv/qBessel.hpp>
#include <cmath>
typedef kv::interval<double> itv;
typedef kv::complex< kv::interval<double> > cp;
using namespace std;
int main()
{
  cout.precision(17);
  int n=20;
  itv x,nu,q,xx;
  q="0.3";
  nu=1.5;
  x=(2.5,2.6);
  xx=mid(x);
  for(int i=1;i<=n;i++){
    xx=xx-kv::Hahn_Exton(itv(xx),itv(nu),itv(q))*(1-q)*x
      /(kv::Hahn_Exton(itv(x),itv(nu),itv(q))-kv::Hahn_Exton(itv(q*x),itv(nu),itv(q)));
  cout<<xx<<endl;
  cout<<"value of HE"<<kv::Hahn_Exton(itv(xx),itv(nu),itv(q))<<endl;
  }
 
}
