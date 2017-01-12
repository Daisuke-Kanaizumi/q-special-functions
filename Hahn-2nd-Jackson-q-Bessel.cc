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
  int K=2500;
  int j=1;
  kv::interval<double> q;
  kv::interval<double> nu; // nu is a real number
  kv::interval<double> x;  
  q=0.1;// q must be positive
  nu=3.5;
  x=1.5;
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
    a = a+j*pow(x/2,2*k)/(b*c);
  }
  b=b*(1-pow(q,K+1));
  c=c*(1-pow(q,K+nu+1));
  double rad=abs(pow(-1,K+1)*pow(x/2,2*K+2)/(b*c)).upper();
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
  kv::interval<double> first=pow(x/2,nu)*series/denom; // first OK
  // from below, Hahn`s formula will be used
  // this formula describes the relation between 1st and 2nd Jackson q-Bessel function
  // calculate ((-x^2)/4;q)_Infinity with q-exponential E_q(z)
  kv::interval<double> z1;
  z1=-pow(x,2)/4;
  kv::interval<double> d;
  kv::interval<double> e;
  d=1/(1-z1);
  e=1.;
  int j1;
  j1=1;
  // derived from the definition of q-Pochhammer symbol
  for (int n=1; n<=K; n++){
    j1 = -1*j1;
    e=e*(1-pow(q,n));
    d=d+j1*pow(q,n*(n+1)/2)/(e*(1-z1*pow(q,n)));
  }
  e=e*(1-pow(q,K+1));
  double rad3=abs(pow(-1,K+1)*pow(q,(K+1)*(K+2)/2)/(e*(1-z1*pow(q,K+1)))).upper();
  kv::interval<double> series2=d+rad3; // result must be an interval
  // series 2 OK
  // calculate the value of Euler function
  kv::interval<double> h=1-q;
  int m1=1;
  for (int r=1; r<=K; r++){
    m1 = -1*m1;
    h = h+m1*(1-pow(q,2*r+1))*pow(q,r*(3*r+1)/2);
  }
  double eulerrad=abs(pow(-1,K+1)*(1-pow(q,2*K+3))*pow(q,(K+1)*(3*K+4)/2)).upper();
  kv::interval<double> euler=h+eulerrad; // result must be an interval
  // euler OK
  kv::interval<double> hahn;
  hahn=euler/series2;
  // hahn OK
  
  kv::interval<double> res;
  res=first*hahn;
  cout << fixed << setprecision(16) <<res << endl;
}
// definition of the q-exponential function is given in the reference below
// Basic Hypergeometric Series Cambridge University Press (Gasper, Rahman, 1990)

// this program needs to be modified
