#include <kv/Heine.hpp>
#include <kv/qAiry.hpp>
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
  int n=20;
  itv nu,q;
  ub::vector< itv > x(100);
  q="0.7";
  nu=1.5;
  x(0)=4.5;
  x(1)=4.45;
  for(int i=1;i<=n;i++){
    x(i+1)=x(i)-kv::Jackson2(itv(x(i)),itv(nu),itv(q))*(x(i)-x(i-1))
      /(kv::Jackson2(itv(x(i)),itv(nu),itv(q))-kv::Jackson2(itv(x(i-1)),itv(nu),itv(q)));
  cout<<x(i+1)<<endl;
  cout<<"value of J2 inf"<<kv::Jackson2(itv(x(i+1).lower()),itv(nu),itv(q))<<endl;
  cout<<"value of J2 sup"<<kv::Jackson2(itv(x(i+1).upper()),itv(nu),itv(q))<<endl;
  cout<<"value of J2 mid"<<kv::Jackson2(itv(mid(x(i+1))),itv(nu),itv(q))<<endl;
  }
}
