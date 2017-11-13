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
  int n=50;
  itv x,nu,q,p;
  p="0.8";
  q="0.3";
  nu=1.5;
  x=2.5;
  for(int i=1;i<=n;i++){
    x=x-kv::Hahn_Exton(itv(x),itv(nu),itv(q))*(p-q)*x
      /(kv::Hahn_Exton(itv(p*x),itv(nu),itv(q))-kv::Hahn_Exton(itv(q*x),itv(nu),itv(q)));
  cout<<x<<endl;
  cout<<"value of HE"<<kv::Hahn_Exton(itv(x),itv(nu),itv(q))<<endl;
  }
 
}
