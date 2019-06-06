#include <kv/Heine.hpp>
#include <kv/qAiry.hpp>
#include <kv/qBessel.hpp>
#include <cmath>
#include <iostream>
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
  itv nu,q,xx;
  q="0.7";
  nu=2.5;
  x(0)=4.5;
  
  for(int i=1;i<=n;i++){
    xx=mid(x(i-1));
    x(i)=xx-kv::Hahn_Exton(itv(xx),itv(nu),itv(q))*(1-q)*x(i-1)
      /(kv::Hahn_Exton(itv(x(i-1)),itv(nu),itv(q))-q*kv::Hahn_Exton(itv(q*x(i-1)),itv(nu),itv(q)));
    
    cout<<x(i)<<endl;
    cout<<"value of HE inf"<<kv::Hahn_Exton(itv(x(i).lower()),itv(nu),itv(q))<<endl;
    cout<<"value of HE sup"<<kv::Hahn_Exton(itv(x(i).upper()),itv(nu),itv(q))<<endl;
    cout<<"value of HE mid"<<kv::Hahn_Exton(itv(mid(x(i))),itv(nu),itv(q))<<endl;

  }
 
}
