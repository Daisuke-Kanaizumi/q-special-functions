#include <kv/Heine.hpp>
#include <kv/qAiry.hpp>
#include <kv/qLaguerre.hpp>
#include <cmath>
typedef kv::interval<double> itv;
typedef kv::complex< kv::interval<double> > cp;
using namespace std;
int main()
{
  cout.precision(17);
  int m=20;
  int n=5;
  itv x,q,alpha;
  q="0.8";
  int nn;
  nn=15;
  x=2.5;
  alpha=3.5;
  for(int i=1;i<=nn;i++){
    x=x-kv::qLpoly(itv(x),itv(q),int(n),itv(alpha))*(1-q)*x
      /(kv::qLpoly(itv(x),itv(q),int(n),itv(alpha))-kv::qLpoly(itv(q*x),itv(q),int(n),itv(alpha)));
  cout<<x<<endl;
  cout<<"value of qLpoly inf"<<kv::qLpoly(itv(x.lower()),itv(q),int(n),itv(alpha))<<endl;
  cout<<"value of qLpoly sup"<<kv::qLpoly(itv(x.upper()),itv(q),int(n),itv(alpha))<<endl;
  cout<<"value of qLpoly mid"<<kv::qLpoly(itv(mid(x)),itv(q),int(n),itv(alpha))<<endl;

  }

}
