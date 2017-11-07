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
  q="0.2";
  nu=1.5;
  x=2.5;
  for(int i=1;i<=n;i++){
    x=x-kv::modified_qBesselI2(itv(x),itv(nu),itv(q))*(1-q)*x
      /(kv::modified_qBesselI2(itv(x),itv(nu),itv(q))-kv::modified_qBesselI2(itv(q*x),itv(nu),itv(q)));
  cout<<x<<endl;
  cout<<"value of I2"<<kv::modified_qBesselI2(itv(x),itv(nu),itv(q))<<endl;
  }
}
