#include <kv/Heine.hpp>
#include <kv/qAiry.hpp>

#include <cmath>
typedef kv::interval<double> itv;
typedef kv::complex< kv::interval<double> > cp;
using namespace std;
int main()
{
  cout.precision(17);
  int n=15;
  itv x,q,xx;
  q="0.9";

  x=2.5;
  for(int i=1;i<=n;i++){
    xx=mid(x);
    x=xx-kv::Ramanujan_qAiry(itv(q),itv(xx))*(1-q)*xx
      /(kv::Ramanujan_qAiry(itv(q),itv(xx))-kv::Ramanujan_qAiry(itv(q),itv(q*xx)))
      +(1-xx*(kv::Ramanujan_qAiry(itv(q),itv(x))-kv::Ramanujan_qAiry(itv(q),itv(q*x))) 
	/x/(kv::Ramanujan_qAiry(itv(q),itv(xx))-kv::Ramanujan_qAiry(itv(q),itv(q*xx))))*(x-xx);
  cout<<x<<endl;
  cout<<"value of RqA inf"<<kv::Ramanujan_qAiry(itv(q),itv(x.lower()))<<endl;
  cout<<"value of RqA sup"<<kv::Ramanujan_qAiry(itv(q),itv(x.upper()))<<endl;
  cout<<"value of RqA mid"<<kv::Ramanujan_qAiry(itv(q),itv(mid(x)))<<endl;
  }
}
