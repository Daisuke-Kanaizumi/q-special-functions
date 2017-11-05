#include <kv/Heine.hpp>
#include <kv/qAiry.hpp>
#include <kv/qBessel.hpp>
#include <cmath>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
typedef kv::interval<double> itv;
typedef kv::complex< kv::interval<double> > cp;
using namespace std;
int main()
{
  cout.precision(17);
  itv x,nu,q,qmvf,xx,qmvf2,qmvf3;
  q="0.3";
  nu=3.5;
  x=itv(5.5,5.6);
  xx=mid(x);
  qmvf=kv::Jackson2(itv(xx),itv(nu),itv(q))
    +(x-xx)*(kv::Jackson2(itv(x),itv(nu),itv(q))-kv::Jackson2(itv(q*x),itv(nu),itv(q)))/(1-q)/x;
  cout<<"value of J2"<<kv::Jackson2(itv(x),itv(nu),itv(q))<<endl;
  cout<<"value of J2 with q-mean value form"<<qmvf<<endl;
  qmvf2=kv::Hahn_Exton(itv(xx),itv(nu),itv(q))
    +(x-xx)*(kv::Hahn_Exton(itv(x),itv(nu),itv(q))-kv::Hahn_Exton(itv(q*x),itv(nu),itv(q)))/(1-q)/x;
  cout<<"value of HE"<<kv::Hahn_Exton(itv(x),itv(nu),itv(q))<<endl;
  cout<<"value of HE with q-mean value form"<<qmvf2<<endl;
  qmvf3=kv::Ramanujan_qAiry(itv(q),itv(xx))
    +(x-xx)*(kv::Ramanujan_qAiry(itv(q),itv(x))-kv::Ramanujan_qAiry(itv(q),itv(q*x)))/(1-q)/x;
  cout<<"value of RqA"<<kv::Ramanujan_qAiry(itv(q),itv(x))<<endl;
  cout<<"value of RqA with q-mean value form"<<qmvf3<<endl;

}
