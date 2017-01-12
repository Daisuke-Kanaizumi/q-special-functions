// verification program for the quantum dilogarithm
// q-extension of dilogarithm
// studied by Faddeev-Kashaev, and Kirillov

// reference
// Don Zagier
// The Dilogarithm Function
// Chapter II, Proposition 2 (power series expansion)

#include<iostream>
#include<cmath>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
using namespace std;
int main(void)
{
cout.precision(16);
  int K=1300;
  int j=1;
  kv::interval<double> q;//-1<q<1
  kv::interval<double> z;//|z|<1 
  q=0.5;
  z=0.1;
  kv::interval<double> a;
  kv::interval<double> b;
  a=1/(1-z);
  b=1.;
 // derived from the definition of q-Pochhammer symbol
  for (int k=1; k<=K; k++){
    j = -1*j;
    b=b*(1-pow(q,k));
    a=a+j*pow(q,k*(k+1)/2)/(b*(1-z*pow(q,k)));
  }
  b=b*(1-pow(q,K+1));
  double rad=abs(pow(-1,K+1)*pow(q,(K+1)*(K+2)/2)/(b*(1-z*pow(q,K+1)))).upper();
  kv::interval<double> series=a+rad; // result must be an interval
  // calculate the value of Euler function
 kv::interval<double> c=1-q;
 int m=1;
 for (int l=1; l<=K; l++){
   m = -1*m;
   c = c+m*(1-pow(q,2*l+1))*pow(q,l*(3*l+1)/2);
  }
  double eulerrad=abs(pow(-1,K+1)*(1-pow(q,2*K+3))*pow(q,(K+1)*(3*K+4)/2)).upper();
  kv::interval<double> euler=c+eulerrad; //final result must be an interval
  kv::interval<double> qp=euler/series;
  kv::interval<double> res=-log(qp);
  cout << res << endl;
}

