#include <kv/Heine.hpp>
#include <kv/qAiry.hpp>
#include <cmath>
typedef kv::interval<double> itv;
typedef kv::complex< kv::interval<double> > cp;
using namespace std;
int main()
{
  cout.precision(17);
  int n=25;
  itv x,nu,q,x0,y;
  q="0.7";
  x=15.;
  x0=x;
  y=1.5;
  for(int i=1;i<=n;i++){
    x=x-(kv::Ramanujan_qAiry(itv(q),itv(x))-y)*(1-q)*x0
      /(kv::Ramanujan_qAiry(itv(q),itv(x0))- q*kv::Ramanujan_qAiry(itv(q),itv(q*x0)));
  cout<<x<<endl;
  cout<<"value of RqA inf"<<kv::Ramanujan_qAiry(itv(q),itv(x.lower()))-y<<endl;
  cout<<"value of RqA sup"<<kv::Ramanujan_qAiry(itv(q),itv(x.upper()))-y<<endl;
  cout<<"value of RqA mid"<<kv::Ramanujan_qAiry(itv(q),itv(mid(x)))-y<<endl;

  }
 
}
