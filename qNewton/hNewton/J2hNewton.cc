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
  itv x,nu,q,h;
  h="0.1";
  q="0.5";
  nu=1.5;
  x=3.3;
  for(int i=1;i<=n;i++){
    x=x-kv::Jackson2(itv(x),itv(nu),itv(q))*h
      /(kv::Jackson2(itv(x+h),itv(nu),itv(q))-kv::Jackson2(itv(x),itv(nu),itv(q)));
  cout<<x<<endl;
  cout<<"value of J2 inf"<<kv::Jackson2(itv(x.lower()),itv(nu),itv(q))<<endl;
  cout<<"value of J2 sup"<<kv::Jackson2(itv(x.upper()),itv(nu),itv(q))<<endl;
  cout<<"value of J2 mid"<<kv::Jackson2(itv(mid(x)),itv(nu),itv(q))<<endl;
  }
  // A. L. Garcia, Numerical Methods for Physics, 2nd Edition, Pearson, 2015.
}
