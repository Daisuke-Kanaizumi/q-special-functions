#include <kv/Heine.hpp>
#include <kv/qAiry.hpp>
#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
typedef kv::interval<double> itv;
typedef kv::complex< kv::interval<double> > cp;
using namespace std;
namespace ub = boost::numeric::ublas;
int main()
{
  cout.precision(17);
  ub::vector< itv > x(100);
  int n=40;
  itv q,xx,xy,omega;
  q="0.9";
  omega="0.3";
  x(0)=1.55;
  xy=mid(x(0));
  for(int i=0;i<=n;i++){
    xx=mid(x(i));
    x(i+1)=xx-kv::Ramanujan_qAiry(itv(q),itv(xx))*(xy-q*xy-omega)
      /(kv::Ramanujan_qAiry(itv(q),itv(xy))-kv::Ramanujan_qAiry(itv(q),itv(q*xy+omega)))
      +(1-(xy-q*xy-omega)*(kv::Ramanujan_qAiry(itv(q),itv(x(i)))-kv::Ramanujan_qAiry(itv(q),itv(q*x(i)+omega))) 
	/(x(i)-q*x(i)-omega)/(kv::Ramanujan_qAiry(itv(q),itv(xy))-kv::Ramanujan_qAiry(itv(q),itv(q*xy+omega))))*(x(i)-xx);
    cout<<x(i+1)<<endl;
    cout<<"value of RqA inf"<<kv::Ramanujan_qAiry(itv(q),itv(x(i+1).lower()))<<endl;
    cout<<"value of RqA sup"<<kv::Ramanujan_qAiry(itv(q),itv(x(i+1).upper()))<<endl;
    cout<<"value of RqA mid"<<kv::Ramanujan_qAiry(itv(q),itv(mid(x(i+1))))<<endl;
    cout<<"value of contraction const:"<<abs(1-(xy-q*xy-omega)*(kv::Ramanujan_qAiry(itv(q),itv(x(i)))-kv::Ramanujan_qAiry(itv(q),itv(q*x(i)+omega))) 
					     /(x(i)-q*x(i)-omega)/(kv::Ramanujan_qAiry(itv(q),itv(xy))-kv::Ramanujan_qAiry(itv(q),itv(q*xy+omega)))).upper()<<endl;
  }
}
