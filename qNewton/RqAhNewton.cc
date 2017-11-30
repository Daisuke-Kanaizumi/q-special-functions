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
  itv x,q,h;
  q="0.9";
  h="0.1";
  x=2.5;
  for(int i=1;i<=n;i++){
    x=x-kv::Ramanujan_qAiry(itv(q),itv(x))*h
      /(kv::Ramanujan_qAiry(itv(q),itv(x+h))-kv::Ramanujan_qAiry(itv(q),itv(x)));
  cout<<x<<endl;
  cout<<"value of RqA inf"<<kv::Ramanujan_qAiry(itv(q),itv(x.lower()))<<endl;
  cout<<"value of RqA sup"<<kv::Ramanujan_qAiry(itv(q),itv(x.upper()))<<endl;
  cout<<"value of RqA mid"<<kv::Ramanujan_qAiry(itv(q),itv(mid(x)))<<endl;
  }
}
