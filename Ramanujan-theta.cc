#include<iostream>
#include<cmath>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
using namespace std;
int main(void)
{
  // verification program for the Ramanujan theta function
  // defined as f(a,b)=Sum[(a^(n*(n+1)/2))*b^(n*(n-1)/2),{n,-Infinity,Infinity}]
  // this function is the generalized form of Jacobi theta functions
  // this function is defined when |ab|<1
  // Jacobi triple product identity is used for verification
  // f(a,b)=(-a,ab)_Infinity (-b,ab)_Infinity (ab,ab)_Infinity
  // f(a,b) is expressed with q exponential functions
cout.precision(16); 
 int K=1000;
 int j=1;
 kv::interval<double> A;
 kv::interval<double> B;
 A=0.1;
 B=0.2;
 // calculate (-a,ab)_Infinity
 kv::interval<double> q1;
 kv::interval<double> z1;
 z1=-A;
 q1=A*B; 
 kv::interval<double> a;
 kv::interval<double> b;
  a=1/(1-z1);
  b=1.;
 // derived from the definition of q-Pochhammer symbol
  for (int k=1; k<=K; k++){
    j = -1*j;
    b=b*(1-pow(q1,k));
    a=a+j*pow(q1,k*(k+1)/2)/(b*(1-z1*pow(q1,k)));
  }
  b=b*(1-pow(q1,K+1));
  double rad1=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b*(1-z1*pow(q1,K+1)))).upper();
  kv::interval<double> series1=a+rad1; // result must be an interval
  // calculate the value of Euler function
 kv::interval<double> c=1-q1;
 int m=1;
 for (int l=1; l<=K; l++){
   m = -1*m;
   c = c+m*(1-pow(q1,2*l+1))*pow(q1,l*(3*l+1)/2);
  }
  double eulerrad1=abs(pow(-1,K+1)*(1-pow(q1,2*K+3))*pow(q1,(K+1)*(3*K+4)/2)).upper();
  kv::interval<double> euler=c+eulerrad1; 
  kv::interval<double> exp1=series1/euler;
  kv::interval<double> pro1=1/exp1;
 // calculate (-b,ab)_Infinity
 kv::interval<double> z2;
 z2=-B;
 // reuse q1
 kv::interval<double> a2;
 kv::interval<double> b2;
  a2=1/(1-z2);
  b2=1.;
 // derived from the definition of q-Pochhammer symbol
 int j2=1;
  for (int k2=1; k2<=K; k2++){
    j2 = -1*j2;
    b2=b2*(1-pow(q1,k2));
    a2=a2+j2*pow(q1,k2*(k2+1)/2)/(b2*(1-z2*pow(q1,k2)));
  }
  b2=b2*(1-pow(q1,K+1));
  double rad2=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b2*(1-z2*pow(q1,K+1)))).upper();
  kv::interval<double> series2=a2+rad2; // result must be an interval
  // reuse the value of Euler function
 
  kv::interval<double> exp2=series2/euler;
  kv::interval<double> pro2=1/exp2;
// reuse (ab,ab)_Infinity (Euler function) 
  kv::interval<double> res=pro1*pro2*euler;
  cout << res << endl;// final result
}
