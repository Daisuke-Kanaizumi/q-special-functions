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
  int n=50;
  itv x,nu,q;
  q="0.7";
  nu=2.5;
  x=3.5;
  for(int i=1;i<=n;i++){
    x=x-kv::Jackson2(itv(x),itv(nu),itv(q))
      /kv::Jackson2(itv(x),itv(nu-1.),itv(q))*(1-q)*x
      /(kv::Jackson2(itv(x),itv(nu),itv(q))/kv::Jackson2(itv(x),itv(nu-1.),itv(q))
	-kv::Jackson2(itv(q*x),itv(nu),itv(q))/kv::Jackson2(itv(q*x),itv(nu-1.),itv(q)));
  cout<<x<<endl;
  cout<<"value of J2"<<kv::Jackson2(itv(x),itv(nu),itv(q))<<endl;
  }
  // Ismail, Mourad E. H. (1982), Journal of Mathematical Analysis and Applications, 86 (1): 1-19
  // Gil, A., Segura, J. and Temme, N.M., 2007. Numerical methods for special functions. Society for Industrial and Applied Mathematics.
}
