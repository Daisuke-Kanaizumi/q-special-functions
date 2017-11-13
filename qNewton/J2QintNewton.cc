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
  itv x,nu,q,xx;
  q="0.3";
  nu=1.5;
  x=itv(2.5,2.6);
  xx=mid(x);
  for(int i=1;i<=n;i++){
    xx=xx-kv::Jackson2(itv(xx),itv(nu),itv(q))*(1-q)*x
      /(kv::Jackson2(itv(x),itv(nu),itv(q))-kv::Jackson2(itv(q*x),itv(nu),itv(q)));
  cout<<xx<<endl;
  cout<<"value of J2"<<kv::Jackson2(itv(xx),itv(nu),itv(q))<<endl;
  }
  // A. L. Garcia, Numerical Methods for Physics, 2nd Edition, Pearson, 2015.
}
