#include <kv/Heine.hpp>
#include <kv/qAiry.hpp>
#include <cmath>
typedef kv::interval<double> itv;
typedef kv::complex< kv::interval<double> > cp;
using namespace std;
int main()
{
  cout.precision(17);
  int n=50;
  itv x,q;
  q="0.7";
 
  x=3.5;
  for(int i=1;i<=n;i++){
    x=x-kv::Ramanujan_qAiry(itv(q),itv(x))*(1-q)*x
      /(kv::Ramanujan_qAiry(itv(q),itv(x))-kv::Ramanujan_qAiry(itv(q),itv(q*x)));
  cout<<x<<endl;
  cout<<"value of RqA inf"<<kv::Ramanujan_qAiry(itv(q),itv(x.lower()))<<endl;
  cout<<"value of RqA sup"<<kv::Ramanujan_qAiry(itv(q),itv(x.upper()))<<endl;
  cout<<"value of RqA mid"<<kv::Ramanujan_qAiry(itv(q),itv(mid(x)))<<endl;
  }
}
