#include<iostream>
#include<cmath>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include<iomanip> // include this to set digits
using namespace std;
int main(void)
{
  // verfication program for the 2nd Jackson  q-Bessel function
  // properties of the q-Pochhammer symbol were used  
  int K=1000;
  int j=1;
  kv::interval<double> q;
  kv::interval<double> nu; // nu is a real number
  kv::interval<double> x;  
  q=0.1;// q must be positive
  nu=3.5;
  x=1.4;
  kv::interval<double> a;
  kv::interval<double> b;
  kv::interval<double> c;
  a=1.;
  b=1.;
  c=1.;
 // derived from the definition of q-Pochhammer symbol
  for (int k=1; k<=K; k++){
    j = -1*j;
    b=b*(1-pow(q,k));
    c=c*(1-pow(q,k+nu));
    a = a+j*pow(q,k*(k-1))*pow(pow(q,nu+1)*x*x/4,k)/(b*c);
  }
  b=b*(1-pow(q,K+1));
  c=c*(1-pow(q,K+nu+1));
  double rad=abs(pow(q,K*(K+1))*pow(pow(q,nu+1)*x*x/4,K+1)/(b*c)).upper();
  // .upper() allows to output supremum
  kv::interval<double> series=a+rad; // result must be an interval
  
  kv::interval<double> z;
  z=pow(q,nu+1);
  kv::interval<double> f;
  kv::interval<double> g;
  f=1/(1-z);
  g=1.;
  int m=1;
 // derived from the definition of q-Pochhammer symbol
  for (int l=1; l<=K; l++){
    m = -1*m;
    g=g*(1-pow(q,l));
    f=f+m*pow(q,l*(l+1)/2)/(g*(1-z*pow(q,l)));
  }
  g=g*(1-pow(q,K+1));
  double rad2=abs(pow(-1,K+1)*pow(q,(K+1)*(K+2)/2)/(g*(1-z*pow(q,K+1)))).upper();
  kv::interval<double> denom=f+rad2; // result must be an interval
  kv::interval<double> res=pow(x/2,nu)*series/denom;
  cout << fixed << setprecision(16) <<res << endl;
}
