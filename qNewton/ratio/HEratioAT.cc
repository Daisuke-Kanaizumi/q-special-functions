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
  itv x,nu,q;
  q="0.7";
  nu=2.5;
  x=1.2;

  for(int i=1;i<=n;i++){
    x=x-atan(kv::HEratio(itv(x),itv(nu),itv(q)));
  cout<<x<<endl;
  cout<<"value of HE inf"<<kv::Hahn_Exton(itv(x.lower()),itv(nu),itv(q))<<endl;
  cout<<"value of HE sup"<<kv::Hahn_Exton(itv(x.upper()),itv(nu),itv(q))<<endl;
  cout<<"value of HE mid"<<kv::Hahn_Exton(itv(mid(x)),itv(nu),itv(q))<<endl;
  }
}
