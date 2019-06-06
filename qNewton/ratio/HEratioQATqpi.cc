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
  q="0.8";
  nu=4.5;
  x=4.5;

  for(int i=1;i<=n;i++){
    x=x-q_atan_qpi(itv(q),itv(kv::HEratio(itv(x),itv(nu),itv(q))));
  cout<<x<<endl;
  cout<<"value of HE inf"<<kv::Hahn_Exton(itv(x.lower()),itv(nu),itv(q))<<endl;
  cout<<"value of HE sup"<<kv::Hahn_Exton(itv(x.upper()),itv(nu),itv(q))<<endl;
  cout<<"value of HE mid"<<kv::Hahn_Exton(itv(mid(x)),itv(nu),itv(q))<<endl;
  }
}
