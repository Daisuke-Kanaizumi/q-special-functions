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
    x(i+1)=x(i)-kv::Hahn_Exton(itv(x(i)),itv(nu),itv(q))*(x(i)-x(i-1))
      /(kv::Hahn_Exton(itv(x(i)),itv(nu),itv(q))-kv::Hahn_Exton(itv(x(i-1)),itv(nu),itv(q)));
  cout<<x(i+1)<<endl;
  cout<<"value of HE"<<kv::Hahn_Exton(itv(x(i+1)),itv(nu),itv(q))<<endl;
  }
}
