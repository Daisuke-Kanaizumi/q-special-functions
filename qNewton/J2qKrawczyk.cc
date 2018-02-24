#include <kv/Heine.hpp>
#include <kv/qBessel.hpp>
#include <cmath>
#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
typedef kv::interval<double> itv;
typedef kv::complex< kv::interval<double> > cp;
using namespace std;
namespace ub = boost::numeric::ublas;
int main()
{
  cout.precision(17);
  ub::vector< itv > x(30);
  int n=20;
  itv nu,q,xx,xy;
  q="0.5";
  nu=1.5;
  x(0)=5.;
  xy=mid(x(0));
  for(int i=0;i<=n;i++){
    xx=mid(x(i));
    x(i+1)=xx-kv::Jackson2(itv(xx),itv(nu),itv(q))*(1-q)*xy
      /(kv::Jackson2(itv(xy),itv(nu),itv(q))-kv::Jackson2(itv(q*xy),itv(nu),itv(q)))
      +(1-xy*(kv::Jackson2(itv(x(i)),itv(nu),itv(q))-kv::Jackson2(itv(q*x(i)),itv(nu),itv(q))) 
	/x(i)/(kv::Jackson2(itv(xy),itv(nu),itv(q))-kv::Jackson2(itv(q*xy),itv(nu),itv(q))))*(x(i)-xx);
    cout<<x(i+1)<<endl;
    cout<<"value of J2 inf"<<kv::Jackson2(itv(x(i+1).lower()),itv(nu),itv(q))<<endl;
    cout<<"value of J2 sup"<<kv::Jackson2(itv(x(i+1).upper()),itv(nu),itv(q))<<endl;
    cout<<"value of J2 mid"<<kv::Jackson2(itv(mid(x(i+1))),itv(nu),itv(q))<<endl;
      cout<<"value of contraction const."<<(abs(1-xy*(kv::Jackson2(itv(x(i+1)),itv(nu),itv(q))-kv::Jackson2(itv(q*x(i+1)),itv(nu),itv(q))) 
					      /x(i)/(kv::Jackson2(itv(xy),itv(nu),itv(q))-kv::Jackson2(itv(q*xy),itv(nu),itv(q))))).upper()<<endl;  
  }
}
