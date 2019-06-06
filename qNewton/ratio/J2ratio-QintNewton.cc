#include <kv/Heine.hpp>
#include <kv/qAiry.hpp>
#include <kv/qBessel.hpp>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/interval.hpp>

// define rounding operations for double
// if "rdouble.hpp" is not included, result of computation is not "verified" 
#include <kv/rdouble.hpp>
typedef kv::interval<double> itv;
typedef kv::complex< kv::interval<double> > cp;
using namespace std;
namespace ub = boost::numeric::ublas;

int main()
{
  cout.precision(17);
  ub::vector< itv > x(30);
  int n=20;
  itv y,nu,q,xx;
  q="0.7";
  nu=1.5;
  y=4.8;
 
  x(0)=y;
  
  for(int i=1;i<=n;i++){
    xx=mid(x(i-1));
    x(i)=xx-kv::J2ratio(itv(xx),itv(nu),itv(q))*(1-q)*x(i-1)
      /(kv::J2ratio(itv(x(i-1)),itv(nu),itv(q))-kv::J2ratio(itv(q*x(i-1)),itv(nu),itv(q)));
    cout<<x(i)<<endl;
    //cout<<i<<endl;
    cout<<"value of J2 inf"<<kv::Jackson2(itv(x(i).lower()),itv(nu),itv(q))<<endl;
    cout<<"value of J2 sup"<<kv::Jackson2(itv(x(i).upper()),itv(nu),itv(q))<<endl;
    cout<<"value of J2 mid"<<kv::Jackson2(itv(mid(x(i))),itv(nu),itv(q))<<endl;
  }
 }
  // to be updated

