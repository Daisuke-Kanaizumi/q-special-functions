#include <kv/Pochhammer.hpp>
#include <kv/Sokal.hpp>
#include <kv/qPochhammerVer2.hpp>
#include <cmath>
#include <kv/Gatteschi.hpp>
// finding fixed point of (z;q) with Steffensen iteration
using namespace std;
typedef kv::interval<double> itv;
typedef kv::complex< kv::interval<double> > cp;
int main()
{
cout.precision(17);
 itv q,z;
 q="0.6";
 z=1.5;
 int n;
 n=20;
 /*for(int i=1;i<=n;i++){
   z=kv::qPVer2(itv(kv::qPVer2(itv(z), itv(q))), itv(q))
     -pow(kv::qPVer2(itv(kv::qPVer2(itv(z), itv(q))), itv(q))-kv::qPVer2(itv(z), itv(q)),2.)
     /(kv::qPVer2(itv(kv::qPVer2(itv(z), itv(q))), itv(q))-2*kv::qPVer2(itv(z), itv(q))+z);
   cout<<z<<endl;
   cout << "value of qp"<<kv::qPVer2(itv(z), itv(q)) << "\n";
   }*/
 for(int i=1;i<=n;i++){
   z=kv::Gatteschi_qp(itv(kv::Gatteschi_qp(itv(z), itv(q))), itv(q))
     -pow(kv::Gatteschi_qp(itv(kv::Gatteschi_qp(itv(z), itv(q))), itv(q))-kv::Gatteschi_qp(itv(z), itv(q)),2.)
     /(kv::Gatteschi_qp(itv(kv::Gatteschi_qp(itv(z), itv(q))), itv(q))-2*kv::Gatteschi_qp(itv(z), itv(q))+z);
   cout<<z<<endl;
   cout << "value of (z;q)-z:"<<kv::qPVer2(itv(z), itv(q))-z << "\n";
 }
}
// q=0.1~0.7,z=1.5
