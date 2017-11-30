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
  ub::vector< itv > x(30);
  int n=25;
  itv nu,q,xx,xy;
  q="0.7";
  nu=2.5;
  x(0)="0.9";
  xy=mid(x(0));
  for(int i=0;i<=n;i++){
    xx=mid(x(i));
    x(i+1)=xx-kv::Hahn_Exton(itv(xx),itv(nu),itv(q))*(1-q)*xy
      /(kv::Hahn_Exton(itv(xy),itv(nu),itv(q))-kv::Hahn_Exton(itv(q*xy),itv(nu),itv(q)))
      +(1-xy*(kv::Hahn_Exton(itv(x(i)),itv(nu),itv(q))-kv::Hahn_Exton(itv(q*x(i)),itv(nu),itv(q))) 
	/x(i)/(kv::Hahn_Exton(itv(xy),itv(nu),itv(q))-kv::Hahn_Exton(itv(q*xy),itv(nu),itv(q))))*(x(i)-xx);
    cout<<x(i+1)<<endl;
    cout<<"value of HE inf"<<kv::Hahn_Exton(itv(x(i+1).lower()),itv(nu),itv(q))<<endl;
    cout<<"value of HE sup"<<kv::Hahn_Exton(itv(x(i+1).upper()),itv(nu),itv(q))<<endl;
    cout<<"value of HE mid"<<kv::Hahn_Exton(itv(mid(x(i+1))),itv(nu),itv(q))<<endl;

  }
}
