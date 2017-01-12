
#include<cmath>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
using namespace std;
int main(void)
{
  // verification program for the Ramanujan q beta integral
  // reference
  // Classical and Quantum Orthogonal Polynomials in One Variable
  // Encyclopedia of Mathematics and its Applications 98
  // Mourad E. H. Ismail
cout.precision(16);
 int K=1000;

 kv::interval<double> A;
 kv::interval<double> B;
 kv::interval<double> q;
 A=0.1;
 B=0.2;
 q=0.3;


   
     //calculate (a,b;q)_Infinity
     // first, (a;q)_Infinity
     int j=1;
     kv::interval<double> q1;
     kv::interval<double> z1;
 z1=A;
 q1=q;
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
  // the value of Euler function will be cancelled

  kv::interval<double> exp1=series1;
  kv::interval<double> pro1=1/exp1;
  // next, calculate  (b;q)_Infinity
  kv::interval<double> z2;
  z2=B;
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
  kv::interval<double> series2=a2+rad2;

  kv::interval<double> exp2=series2;
  kv::interval<double> pro2=1/exp2;
  // calculate (q,qab;q)_Infinity
  // (q,q)_Infinity will be cancelled
  // (ab/q;q)_Infinity
 kv::interval<double> z3;
 z3=A*B/q;
 // reuse q1
 kv::interval<double> a3;
 kv::interval<double> b3;
  a3=1/(1-z3);
  b3=1.;
 // derived from the definition of q-Pochhammer symbol
 int j3=1;
  for (int k3=1; k3<=K; k3++){
    j3 = -1*j3;
    b3=b3*(1-pow(q1,k3));
    a3=a3+j3*pow(q1,k3*(k3+1)/2)/(b3*(1-z3*pow(q1,k3)));
  }
  b3=b3*(1-pow(q1,K+1));
  double rad3=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b3*(1-z3*pow(q1,K+1)))).upper();
  kv::interval<double> series3=a3+rad3; // result must be an interval

  kv::interval<double> exp3=series3;
  kv::interval<double> pro3=1/exp3;

  kv::interval<double> res=(pro1*pro2)/(pro3);
  cout << res << endl;// final result

   }

