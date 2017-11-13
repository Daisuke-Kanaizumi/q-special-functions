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
  itv x,nu,q,p;
  p="0.5";
  q="0.8";
  nu=1.5;
  x=2.5;
  for(int i=1;i<=n;i++){
    x=x-kv::Jackson2(itv(x),itv(nu),itv(q))*(p-q)*x
      /(kv::Jackson2(itv(p*x),itv(nu),itv(q))-kv::Jackson2(itv(q*x),itv(nu),itv(q)));
  cout<<x<<endl;
  cout<<"value of J2"<<kv::Jackson2(itv(x),itv(nu),itv(q))<<endl;
  }
  // A. L. Garcia, Numerical Methods for Physics, 2nd Edition, Pearson, 2015.
}
